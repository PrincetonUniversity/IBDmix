#pragma once
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "IBDmix/Mask_Reader.h"
#include "IBDmix/Sample_Mapper.h"

int find_token(const char * query, const char * str);
const unsigned char IN_MASK = 1 << 0;
const unsigned char MAF_LOW = 1 << 1;
const unsigned char MAF_HIGH = 1 << 2;
const unsigned char RECOVER_2_0 = 1 << 3;
const unsigned char RECOVER_0_2 = 1 << 4;


class Genotype_Reader{
    private:
        FILE * genotype;
        size_t buf_size;
        Mask_Reader mask;
        Sample_Mapper sample_mapper;

    public:
        Genotype_Reader(FILE * genotype, std::istream *mask=nullptr,
                double archaic_error=0.01, double modern_error_max=0.002,
                double modern_error_proportion=2, double minesp=1e-200,
                int minor_allele_cutoff=1);
        ~Genotype_Reader();

        char *buffer;
        int minor_allele_cutoff;
        double archaic_error, modern_error_max, modern_error_proportion,
               minesp;

        std::vector<double> lod_scores, lod_cache;
        std::vector<unsigned char> recover_type;

        void process_line_buffer(bool selected);
        double get_modern_error(double frequency);
        bool get_frequency(double &frequency);

        void update_lod_cache(char archaic, double freq_b, double modern_error,
                bool selected=true);
        double calculate_lod(char modern);

        int chromosome;
        unsigned long int position;
        unsigned char line_filtering;
        int num_samples() { return sample_mapper.size(); }

        int initialize(std::istream &samples, std::string archaic="");
        const std::vector<std::string>& get_samples();
        bool update(void);
};
// TODO make a LOD Calculator class, needs mapper
// TODO handle string chromosomes, deal with sorting on mask
// std::istringstream iss(buffer);
// std::copy(std::istream_iterator<std::string>(iss),
//         std::istream_iterator<std::string>(),
//         std::back_inserter(sample_list));
// result = sample_list.size() - 1;  // remove 1 for archaic
