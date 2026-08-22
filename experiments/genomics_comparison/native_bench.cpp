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
  if (argc != 6) {
    std::cerr << "usage: caramel_bench TABLE none|bloom BPE HASHES REPEATS\n";
    return 2;
  }
  const std::string table = argv[1];
  const std::string kind = argv[2];
  const size_t bpe = std::stoull(argv[3]);
  const size_t hashes = std::stoull(argv[4]);
  const int repeats = std::stoi(argv[5]);
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
  if (kind == "bloom") config = std::make_shared<caramel::BloomPreFilterConfig>(bpe, hashes);
  else if (kind != "none") return 4;

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
    std::mt19937_64 generator(42 + repetition);
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
  std::cout << "{\"serialized_bytes\":" << static_cast<uint64_t>(median(sizes))
            << ",\"build_seconds\":" << median(builds)
            << ",\"query_ns\":" << median(queries)
            << ",\"repetitions\":" << repeats << "}\n";
}
