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
  unsigned char getLineFilter() const { return line_filtering; }
  unsigned char getRecoverType(int index) const { return recover_type[index]; }
  double getLodScore(int index) const { return lod_scores[index]; }
  char getArchaic() const { return archaic; }
  char getAlt() const { return alt; }
  char getRef() const { return ref; }
  uint64_t getPosition() const { return position; }
  double getAlleleFrequency() const { return allele_frequency; }
  const std::string &getChromosome() const { return chromosome; }

 private:
  Mask_Reader mask;
  Sample_Mapper sample_mapper;
  LodCalculator calculator;

  std::istream *genotype;
  std::istringstream iss;
  std::string token;
  std::string buffer;
  std::string chromosome;
  std::vector<unsigned char> recover_type;
  std::vector<double> lod_scores;

  int minor_allele_cutoff;
  unsigned char line_filtering;
  char archaic;
  char alt;
  char ref;
  uint64_t position;
  double allele_frequency = 0;

  bool find_frequency();
  void process_line_buffer(bool selected);
};
