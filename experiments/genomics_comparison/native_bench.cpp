#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "src/construct/Construct.h"
#include "src/construct/filter/FilterConfig.h"

using Clock = std::chrono::steady_clock;

double median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  const auto middle = values.size() / 2;
  return values.size() % 2 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
}

int main(int argc, char** argv) {
  if (argc != 6 && argc != 7) {
    std::cerr << "usage: caramel_bench TABLE none|bloom BPE HASHES BUILDS [QUERY_BATCHES]\n"
                 "       caramel_bench TABLE compare BPE HASHES ROUNDS BATCH_SIZE\n";
    return 2;
  }
  const std::string table = argv[1];
  const std::string kind = argv[2];
  const size_t bpe = std::stoull(argv[3]);
  const size_t hashes = std::stoull(argv[4]);
  const int repeats = std::stoi(argv[5]);
  // Query batches repeat within a build: construction is slow, lookups are noisy.
  const int query_batches = argc == 7 ? std::stoi(argv[6]) : 1;
  std::ifstream input(table);
  if (!input) {
    std::cerr << "cannot open " << table << '\n';
    return 1;
  }
  std::vector<std::string> keys;
  std::vector<uint32_t> values;
  std::string line;
  while (std::getline(input, line)) {
    const auto tab = line.find('\t');
    if (tab == std::string::npos) return 3;
    keys.emplace_back(line, 0, tab);
    values.push_back(static_cast<uint32_t>(std::stoul(line.substr(tab + 1))));
  }
  caramel::PreFilterConfigPtr config = nullptr;
  if (kind == "bloom" || kind == "compare")
    config = std::make_shared<caramel::BloomPreFilterConfig>(bpe, hashes);
  else if (kind != "none") return 4;

  // Both indexes are built in one process and answer the same keys in the same
  // round, alternating which runs first, so machine drift hits both equally.
  if (kind == "compare") {
    const int rounds = std::stoi(argv[5]);
    const size_t batch = std::stoull(argv[6]);
    auto plain = caramel::constructCsf<uint32_t>(keys, values, nullptr, false);
    auto filtered = caramel::constructCsf<uint32_t>(keys, values, config, false);

    auto timed = [&](const auto& structure, const std::vector<size_t>& indexes) {
      volatile uint64_t sink = 0;
      const auto start = Clock::now();
      for (size_t index : indexes) sink += structure->query(keys[index]);
      const auto elapsed = std::chrono::duration<double, std::nano>(Clock::now() - start).count();
      if (sink == 0xdeadbeef) std::cerr << sink;
      return elapsed / indexes.size();
    };

    // Latency, not throughput: chaining each key address to the previous
    // result stops the lookups overlapping. The clock (41.67 ns on Apple
    // Silicon) is far too coarse to time one query directly.
    auto single = [&](const auto& structure, const std::vector<size_t>& indexes) {
      size_t at = 0;
      volatile uint64_t sink = 0;
      const auto start_time = Clock::now();
      for (size_t i = 0; i < indexes.size(); ++i) {
        const uint32_t value = structure->query(keys[indexes[at]]);
        sink += value;
        at = (at + 1 + (value & 7)) % indexes.size();
      }
      const auto elapsed = std::chrono::duration<double, std::nano>(Clock::now() - start_time).count();
      if (sink == 0xdeadbeef) std::cerr << sink;
      return elapsed / indexes.size();
    };

    auto measure = [&](const std::string& path, const auto& structure) {
      structure->save(path, 1);
      const auto bytes = static_cast<double>(std::filesystem::file_size(path));
      std::filesystem::remove(path);
      return bytes;
    };
    const double plain_bytes = measure(table + ".plain.csf", plain);
    const double filter_bytes = measure(table + ".filter.csf", filtered);

    std::vector<double> plain_ns, filter_ns, deltas, plain_single, filter_single;
    int filter_faster = 0;
    for (int round = 0; round < rounds; ++round) {
      // Uniform over the key set, so a fraction alpha of queries carry the
      // dominating value.
      std::mt19937_64 generator(42 + round);
      std::uniform_int_distribution<size_t> sample(0, keys.size() - 1);
      std::vector<size_t> indexes(batch);
      for (auto& index : indexes) index = sample(generator);
      // Warm the whole sample on both structures: warming only a prefix leaves
      // whichever arm runs first paying to fault in the rest.
      volatile uint64_t warm = 0;
      for (size_t index : indexes) warm += plain->query(keys[index]) + filtered->query(keys[index]);
      (void)warm;
      double p, f, ps, fs;
      if (round % 2 == 0) {
        p = timed(plain, indexes); f = timed(filtered, indexes);
        ps = single(plain, indexes); fs = single(filtered, indexes);
      } else {
        f = timed(filtered, indexes); p = timed(plain, indexes);
        fs = single(filtered, indexes); ps = single(plain, indexes);
      }
      plain_single.push_back(ps);
      filter_single.push_back(fs);
      plain_ns.push_back(p);
      filter_ns.push_back(f);
      deltas.push_back(p - f);
      if (f < p) ++filter_faster;
    }
    std::cout << "{\"plain_serialized_bytes\":" << static_cast<uint64_t>(plain_bytes)
              << ",\"filter_serialized_bytes\":" << static_cast<uint64_t>(filter_bytes)
              << ",\"plain_query_ns\":" << median(plain_ns)
              << ",\"filter_query_ns\":" << median(filter_ns)
              << ",\"plain_single_ns\":" << median(plain_single)
              << ",\"filter_single_ns\":" << median(filter_single)
              << ",\"paired_delta_ns\":" << median(deltas)
              << ",\"filter_faster_rounds\":" << filter_faster
              << ",\"rounds\":" << rounds
              << ",\"batch\":" << batch << "}\n";
    return 0;
  }

  std::vector<double> builds, queries, sizes;
  for (int repetition = 0; repetition < repeats; ++repetition) {
    const auto build_start = Clock::now();
    auto structure = caramel::constructCsf<uint32_t>(keys, values, config, false);
    builds.push_back(std::chrono::duration<double>(Clock::now() - build_start).count());
    const std::string serialized = table + ".bench.csf";
    structure->save(serialized, 1);
    sizes.push_back(static_cast<double>(std::filesystem::file_size(serialized)));
    std::filesystem::remove(serialized);
    for (size_t index = 0; index < std::min<size_t>(1000, keys.size()); ++index) {
      if (structure->query(keys[index]) != values[index]) {
        std::cerr << "incorrect query at row " << index << '\n';
        return 5;
      }
    }
    for (int batch = 0; batch < query_batches; ++batch) {
      std::mt19937_64 generator(42 + repetition * query_batches + batch);
      std::uniform_int_distribution<size_t> sample(0, keys.size() - 1);
      std::vector<size_t> indexes(10000);
      for (auto& index : indexes) index = sample(generator);
      volatile uint64_t sink = 0;
      for (size_t index = 0; index < 100; ++index) sink += structure->query(keys[indexes[index]]);
      const auto query_start = Clock::now();
      for (size_t index : indexes) sink += structure->query(keys[index]);
      const auto elapsed = std::chrono::duration<double, std::nano>(Clock::now() - query_start).count();
      queries.push_back(elapsed / indexes.size());
      if (sink == 0xdeadbeef) std::cerr << sink;
    }
  }
  std::cout << "{\"serialized_bytes\":" << static_cast<uint64_t>(median(sizes))
            << ",\"build_seconds\":" << median(builds)
            << ",\"query_ns\":" << median(queries)
            << ",\"query_ns_min\":" << *std::min_element(queries.begin(), queries.end())
            << ",\"query_ns_max\":" << *std::max_element(queries.begin(), queries.end())
            << ",\"repetitions\":" << repeats
            << ",\"query_batches\":" << query_batches << "}\n";
}
