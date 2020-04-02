#pragma once
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "IBDmix/IBD_Segment.h"
#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Genotype_Reader.h"

class IBD_Collection{
    private:
        std::vector<IBD_Segment> IBDs;
        int num_samples;
        IBD_Pool pool;

    public:
        IBD_Collection();
        void initialize(int num_samples, double threshold,
                Genotype_Reader &reader, bool exclusive_end=true,
                bool more_stats=false);
        void update(Genotype_Reader &reader, std::ostream &output);
        void purge(std::ostream &output);
};
