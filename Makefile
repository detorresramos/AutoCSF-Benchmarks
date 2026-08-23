PYTHON := .venv/bin/python

.PHONY: datasets bound-validation synthetic-comparison genomics-comparison \
        reproduce verify baselines publish-datasets

# Genomics inputs: ~500 MB processed, ~10 GB working space. Not needed by the
# two synthetic experiments.
datasets:
	$(PYTHON) datasets/manage.py download
	$(PYTHON) datasets/manage.py validate

bound-validation:
	$(PYTHON) experiments/bound_validation/run.py
	$(PYTHON) experiments/bound_validation/plot.py

synthetic-comparison:
	$(PYTHON) experiments/synthetic_comparison/run.py
	$(PYTHON) experiments/synthetic_comparison/plot.py

genomics-comparison:
	./experiments/genomics_comparison/run.sh
	$(PYTHON) experiments/genomics_comparison/latency.py --dataset rice

reproduce: datasets bound-validation synthetic-comparison genomics-comparison verify

# Compare everything in results/ against the accepted numbers in baselines/ and
# write results/reproduction/receipt.json.
verify:
	$(PYTHON) datasets/manage.py validate
	$(PYTHON) scripts/verify.py

# Accept the numbers currently in results/ as the new baselines.
baselines:
	$(PYTHON) scripts/export_baseline.py

# Stage the genomics tables into an upload folder and push them to the
# Hugging Face dataset repo named in datasets/sources.json.
publish-datasets:
	$(PYTHON) datasets/manage.py stage
	$(PYTHON) datasets/manage.py publish
