#pragma once

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class Sample_Mapper {
 public:
  // read in lines and calls find archaic and map
  // set samples to a nullptr istream if empty
  int initialize(std::istream &genotype, std::istream &samples,
                 std::string archaic = "");

  int size() const { return samples.size(); }
  int getArchaicIndex() const { return archaic_index; }
  int getSample(int index) const { return sample_to_index[index]; }
  const std::vector<std::string> &getSamples() const { return samples; }

 private:
  int archaic_index;
  std::vector<int> sample_to_index;
  std::vector<std::string> samples;

  void map(std::vector<std::string> requested_samples);
  void find_archaic(const std::string &archaic);
};
