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
                Genotype_Reader reader);
        // TODO change output to file after unit testing is "done"
        void update(Genotype_Reader reader, FILE * output);
};
