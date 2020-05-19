#pragma once

#include <string>
#include <vector>

#include "IBDmix/Mask_Reader.h"
#include "IBDmix/Sample_Mapper.h"
#include "IBDmix/lod_calculator.h"

constexpr unsigned char IN_MASK = 1 << 0;
constexpr unsigned char MAF_LOW = 1 << 1;
constexpr unsigned char MAF_HIGH = 1 << 2;
constexpr unsigned char RECOVER_2_0 = 1 << 3;
constexpr unsigned char RECOVER_0_2 = 1 << 4;

class Genotype_Reader {
 public:
  double allele_frequency = 0;

  std::vector<double> lod_scores;
  std::vector<unsigned char> recover_type;

  std::string chromosome;
  uint64_t position;
  char archaic, ref, alt;
  unsigned char line_filtering;

  Genotype_Reader(std::istream *genotype, std::istream *mask = nullptr,
                  double archaic_error = 0.01, double modern_error_max = 0.002,
                  double modern_error_proportion = 2, double minesp = 1e-200,
                  int minor_allele_cutoff = 1)
      : genotype(genotype),
        mask(mask),
        calculator(archaic_error, modern_error_max, modern_error_proportion,
                   minesp),
        minor_allele_cutoff(minor_allele_cutoff) {}

  int initialize(std::istream &samples, std::string archaic = "");
  bool update(void);

  const std::vector<std::string> &get_samples() const;
  int num_samples() const { return sample_mapper.size(); }
  double calculate_lod(char modern) const {
    return calculator.calculate_lod(modern);
  }

 private:
  std::istream *genotype;
  std::istringstream iss;
  std::string token;
  std::string buffer;
  int minor_allele_cutoff;

  Mask_Reader mask;
  Sample_Mapper sample_mapper;
  LodCalculator calculator;

  bool find_frequency();
  void process_line_buffer(bool selected);
};
