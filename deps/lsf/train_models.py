"""Train TFLite models on .lrbin data for use with ribbon_learned_bench.

Reads {dataset}_X.lrbin and {dataset}_y.lrbin, trains MLP models,
exports to TFLite, and writes _eval.txt metadata files.

Usage:
    python train_models.py /data/ caramel
"""

import os
import struct
import sys
import time

import numpy as np


class GlobalMinEarlyStopping:
    """Early stopping from the LSF paper (train/train.py lines 14-55).

    Stops training if the global minimum loss hasn't improved by at least
    `min_improvement_rel` relative to the minimum excluding the last `window` epochs.
    """
    def __init__(self, monitor, min_improvement_rel, window):
        self.monitor = monitor
        self.min_improvement_rel = min_improvement_rel
        self.window = window
        self.losses = []
        self.model = None

    def set_model(self, model):
        self.model = model

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        current_loss = logs.get(self.monitor)
        if current_loss is None:
            return
        self.losses.append(current_loss)
        if len(self.losses) > self.window:
            global_min = min(self.losses)
            min_excl_recent = min(self.losses[:-self.window])
            threshold = min_excl_recent * (1 - self.min_improvement_rel)
            if global_min >= threshold:
                print(
                    f"Stopping at epoch {epoch+1}: global min {global_min:.5f} "
                    f"not {self.min_improvement_rel*100:.0f}% better than "
                    f"min before last {self.window} epochs ({min_excl_recent:.5f})"
                )
                self.model.stop_training = True


def _make_early_stopping_callback(monitor='val_loss'):
    """Create a Keras callback wrapping GlobalMinEarlyStopping."""
    import keras

    class _Callback(keras.callbacks.Callback):
        def __init__(self):
            super().__init__()
            self._es = GlobalMinEarlyStopping(monitor, 0.01, 3)

        def on_epoch_end(self, epoch, logs=None):
            self._es.set_model(self.model)
            self._es.on_epoch_end(epoch, logs)

    return _Callback()


def read_lrbin(data_dir, dataset_name):
    """Read .lrbin feature and label files."""
    x_path = os.path.join(data_dir, f"{dataset_name}_X.lrbin")
    y_path = os.path.join(data_dir, f"{dataset_name}_y.lrbin")

    with open(x_path, "rb") as f:
        num_examples = struct.unpack("<Q", f.read(8))[0]
        num_features = struct.unpack("<Q", f.read(8))[0]
        X = np.frombuffer(f.read(), dtype=np.float32).reshape(num_examples, num_features)

    with open(y_path, "rb") as f:
        num_classes = struct.unpack("<H", f.read(2))[0]
        y = np.frombuffer(f.read(), dtype=np.uint16).astype(np.int32)

    return X, y, num_classes, num_features


def write_lrbin_x(data_dir, dataset_name, X):
    """Overwrite X.lrbin with standardized features."""
    x_path = os.path.join(data_dir, f"{dataset_name}_X.lrbin")
    num_examples, num_features = X.shape
    with open(x_path, "wb") as f:
        f.write(struct.pack("<Q", num_examples))
        f.write(struct.pack("<Q", num_features))
        f.write(X.astype(np.float32).tobytes())


def train_and_export(data_dir, dataset_name, X_train, X_test, y_train, y_test,
                     num_classes, num_layers, hidden_units):
    """Train a single MLP and export to TFLite."""
    import tensorflow as tf
    import keras

    if os.environ.get("LSF_FORCE_CPU"):
        tf.config.set_visible_devices([], "GPU")
    keras.utils.set_random_seed(42)

    model_name = f"{dataset_name}_mlp_L{num_layers}_H{hidden_units}"
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    filename = f"{model_name}_{timestamp}"

    if num_classes > 2:
        output_activation = "softmax"
        output_units = num_classes
        loss = "sparse_categorical_crossentropy"
        metrics = [keras.metrics.SparseCategoricalAccuracy()]
    else:
        output_activation = "sigmoid"
        output_units = 1
        loss = "binary_crossentropy"
        metrics = [keras.metrics.BinaryAccuracy()]

    layers = [keras.layers.Input(shape=(X_train.shape[1],), name="input")]
    for _ in range(num_layers):
        layers.append(keras.layers.Dense(hidden_units, activation="relu"))
    layers.append(keras.layers.Dense(output_units, activation=output_activation))
    model = keras.models.Sequential(layers)

    model.compile(
        loss=loss,
        optimizer=keras.optimizers.Adam(learning_rate=1e-3),
        metrics=metrics,
    )

    print(f"Training {filename}...")
    model.summary()

    early_stop = _make_early_stopping_callback('val_loss')

    t0 = time.perf_counter()
    model.fit(
        X_train, y_train,
        epochs=100,
        batch_size=128,
        verbose=2,
        validation_split=0.1,
        callbacks=[early_stop],
    )
    training_seconds = time.perf_counter() - t0
    model_params = model.count_params()

    model_dir = os.path.join(data_dir, f"{dataset_name}_models")
    os.makedirs(model_dir, exist_ok=True)

    keras_path = os.path.join(model_dir, f"{filename}.keras")
    model.save(keras_path)

    def representative_dataset():
        for x in X_train[:1000]:
            yield {"input": x.reshape(1, -1)}

    # Export TFLite — float16 only (canonical quantization)
    for quantization in ["float16"]:
        converter = tf.lite.TFLiteConverter.from_keras_model(model)
        converter.target_spec.supported_ops = [tf.lite.OpsSet.TFLITE_BUILTINS]
        converter.optimizations = [tf.lite.Optimize.DEFAULT]
        converter.inference_input_type = tf.float32
        converter.inference_output_type = tf.float32
        if quantization == "uint8":
            converter.target_spec.supported_ops = [tf.lite.OpsSet.TFLITE_BUILTINS_INT8]
            converter.representative_dataset = representative_dataset
        elif quantization == "float16":
            converter.target_spec.supported_types = [tf.float16]

        tflite_model = converter.convert()
        tflite_path = os.path.join(model_dir, f"{filename}_{quantization}.tflite")
        with open(tflite_path, "wb") as f:
            f.write(tflite_model)

        # Evaluate (batched to avoid OOM on large datasets)
        interp = tf.lite.Interpreter(model_content=tflite_model)
        interp.allocate_tensors()
        runner = interp.get_signature_runner("serving_default")
        eval_batch = 100_000
        correct = 0
        for start in range(0, len(X_test), eval_batch):
            end = min(start + eval_batch, len(X_test))
            batch_y = runner(input=X_test[start:end])
            batch_y = list(batch_y.values())[0]
            if num_classes == 2:
                batch_y = np.column_stack([1 - batch_y, batch_y])
            correct += np.sum(np.argmax(batch_y, axis=1) == y_test[start:end])
        accuracy = correct / len(X_test)

        eval_path = tflite_path + "_eval.txt"
        with open(eval_path, "w") as f:
            f.write(f"model_l={num_layers} ")
            f.write(f"model_h={hidden_units} ")
            f.write(f"quant={quantization} ")
            f.write(f"training_seconds={training_seconds} ")
            f.write(f"model_params={model_params} ")
            f.write(f"test_accuracy={accuracy * 100}")

        print(f"  Exported {tflite_path} (accuracy={accuracy:.4f})")

    return filename


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <data_dir> <dataset_name>")
        sys.exit(1)

    data_dir = sys.argv[1]
    dataset_name = sys.argv[2]

    print(f"Reading {dataset_name} from {data_dir}...")
    X, y, num_classes, num_features = read_lrbin(data_dir, dataset_name)
    print(f"  {len(X)} examples, {num_features} features, {num_classes} classes")

    # Standardize features (zero mean, unit variance) for better training
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler()
    X = sc.fit_transform(X).astype(np.float32)

    # Overwrite .lrbin with standardized features so ribbon_learned_bench
    # uses the same features the model was trained on
    write_lrbin_x(data_dir, dataset_name, X)

    # Train/test split (stratify only if all classes have >= 2 members)
    from sklearn.model_selection import train_test_split
    from collections import Counter
    min_count = min(Counter(y).values())
    stratify = y if min_count >= 2 else None
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=stratify
    )

    for num_layers, hidden_units in [(0, 0), (1, 50), (1, 100), (2, 50)]:
        print(f"\n{'='*60}")
        print(f"Architecture: L={num_layers}, H={hidden_units}")
        print(f"{'='*60}")
        train_and_export(
            data_dir, dataset_name,
            X_train, X_test, y_train, y_test,
            num_classes, num_layers, hidden_units,
        )

    print("\nTraining complete.")


if __name__ == "__main__":
    main()
