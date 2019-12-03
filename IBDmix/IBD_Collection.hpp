#pragma once
#include <vector>
#include <stdio.h>
#include "IBD_Segment.hpp"
#include "Genotype_Reader.hpp"

class IBD_Collection{
    private:
        std::vector<IBD_Segment> IBDs;
        int num_samples;

    public:
        IBD_Collection();
        void initialize(int num_samples, double threshold,
                Genotype_Reader &reader);
        void update(Genotype_Reader &reader, FILE * output, FILE * sample);
        void purge(FILE * output);
};
