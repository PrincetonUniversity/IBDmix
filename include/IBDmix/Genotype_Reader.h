#pragma once

#include <string>
#include <vector>

#include "IBDmix/Mask_Reader.h"
#include "IBDmix/Sample_Mapper.h"

const unsigned char IN_MASK = 1 << 0;
const unsigned char MAF_LOW = 1 << 1;
const unsigned char MAF_HIGH = 1 << 2;
const unsigned char RECOVER_2_0 = 1 << 3;
const unsigned char RECOVER_0_2 = 1 << 4;

class Genotype_Reader {
 private:
  std::istream *genotype;
  std::istringstream iss;
  std::string token;
  Mask_Reader mask;
  Sample_Mapper sample_mapper;

 public:
  Genotype_Reader(std::istream *genotype, std::istream *mask = nullptr,
                  double archaic_error = 0.01, double modern_error_max = 0.002,
                  double modern_error_proportion = 2, double minesp = 1e-200,
                  int minor_allele_cutoff = 1);

  std::string buffer;
  double archaic_error, modern_error_max, modern_error_proportion, minesp,
      allele_frequency = 0;
  int minor_allele_cutoff;

  std::vector<double> lod_scores, lod_cache;
  std::vector<unsigned char> recover_type;

  void process_line_buffer(bool selected);
  bool find_frequency();
  double get_modern_error();

  void update_lod_cache(char archaic, double freq_b, double modern_error,
                        bool selected = true);
  double calculate_lod(char modern);

  int chromosome;
  uint64_t position;
  char archaic, ref, alt;
  unsigned char line_filtering;
  int num_samples() { return sample_mapper.size(); }

  int initialize(std::istream &samples, std::string archaic = "");
  const std::vector<std::string> &get_samples() const;
  bool update(void);
};
// TODO(troycomi) make a LOD Calculator class, needs mapper
// TODO(troycomi) handle string chromosomes, deal with sorting on mask
