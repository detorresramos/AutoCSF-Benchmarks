#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

namespace {
constexpr unsigned K = 15;
constexpr uint64_t TABLE_SIZE = uint64_t{1} << (2 * K);
constexpr uint64_t MASK = TABLE_SIZE - 1;

int base_code(char value) {
  switch (value) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
  }
}

void consume_sequence(const std::string& sequence, std::vector<uint32_t>& counts) {
  uint64_t word = 0;
  unsigned valid = 0;
  for (char character : sequence) {
    const int code = base_code(character);
    if (code < 0) {
      word = 0;
      valid = 0;
      continue;
    }
    word = ((word << 2) | static_cast<unsigned>(code)) & MASK;
    if (++valid >= K) ++counts[word];
  }
}

bool gz_line(gzFile file, std::string& output) {
  output.clear();
  std::array<char, 1 << 16> buffer{};
  while (gzgets(file, buffer.data(), static_cast<int>(buffer.size())) != nullptr) {
    output += buffer.data();
    if (!output.empty() && output.back() == '\n') break;
  }
  if (output.empty()) return false;
  while (!output.empty() && (output.back() == '\n' || output.back() == '\r')) output.pop_back();
  return true;
}

void count_fasta(gzFile file, std::vector<uint32_t>& counts) {
  std::string line, sequence;
  while (gz_line(file, line)) {
    if (!line.empty() && line.front() == '>') {
      if (!sequence.empty()) consume_sequence(sequence, counts);
      sequence.clear();
    } else {
      sequence += line;
    }
  }
  if (!sequence.empty()) consume_sequence(sequence, counts);
}

void count_fastq(gzFile file, std::vector<uint32_t>& counts) {
  std::string header, sequence, plus, qualities;
  while (gz_line(file, header)) {
    if (!gz_line(file, sequence) || !gz_line(file, plus) || !gz_line(file, qualities)) {
      throw std::runtime_error("truncated FASTQ record");
    }
    if (header.empty() || header.front() != '@' || plus.empty() || plus.front() != '+') {
      throw std::runtime_error("invalid four-line FASTQ record");
    }
    consume_sequence(sequence, counts);
  }
}

std::string decode(uint64_t value) {
  static constexpr char alphabet[] = "ACGT";
  std::string result(K, 'A');
  for (int position = K - 1; position >= 0; --position) {
    result[position] = alphabet[value & 3];
    value >>= 2;
  }
  return result;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: count_kmers STATS_PATH fasta|fastq INPUT...\n";
    return 2;
  }
  try {
    const std::string stats_path = argv[1];
    const std::string format = argv[2];
    auto counts = std::make_unique<std::vector<uint32_t>>(TABLE_SIZE, 0);
    for (int index = 3; index < argc; ++index) {
      gzFile input = gzopen(argv[index], "rb");
      if (!input) throw std::runtime_error("cannot open " + std::string(argv[index]));
      if (format == "fasta") count_fasta(input, *counts);
      else if (format == "fastq") count_fastq(input, *counts);
      else throw std::runtime_error("format must be fasta or fastq");
      gzclose(input);
    }
    std::unordered_map<uint32_t, uint64_t> histogram;
    uint64_t records = 0;
    for (uint64_t word = 0; word < TABLE_SIZE; ++word) {
      const uint32_t count = (*counts)[word];
      if (!count) continue;
      std::cout << decode(word) << '\t' << count << '\n';
      ++histogram[count];
      ++records;
    }
    uint32_t modal_value = 0;
    uint64_t modal_records = 0;
    for (const auto& [value, frequency] : histogram) {
      if (frequency > modal_records || (frequency == modal_records && value < modal_value)) {
        modal_value = value;
        modal_records = frequency;
      }
    }
    std::ofstream stats(stats_path);
    stats << "records\t" << records << '\n'
          << "distinct_values\t" << histogram.size() << '\n'
          << "modal_value\t" << modal_value << '\n'
          << "modal_records\t" << modal_records << '\n';
  } catch (const std::exception& error) {
    std::cerr << "count_kmers: " << error.what() << '\n';
    return 1;
  }
}
