#ifndef INCLUDE_IBDMIX_SAMPLE_MAPPER_H_
#define INCLUDE_IBDMIX_SAMPLE_MAPPER_H_

#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>

class Sample_Mapper{
 public:
    std::vector<int> sample_to_index;
    std::vector<std::string> samples;
    int archaic_index;

    // read in lines and calls find archaic and map
    // set samples to a nullptr istream if empty
    int initialize(std::istream &genotype, std::istream &samples,
            std::string archaic = "");
    // find index of the archaic
    void find_archaic(const std::string &archaic);
    // populate the sample_to_index map
    void map(std::vector<std::string> requested_samples);
    int size() { return samples.size(); }
};

#endif  // INCLUDE_IBDMIX_SAMPLE_MAPPER_H_
