#pragma once
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <memory>
#include "IBDmix/IBD_Segment.h"
#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Genotype_Reader.h"

class IBD_Collection{
    private:
        std::vector<std::unique_ptr<IBD_Segment>> IBDs;
        double threshold;
        bool exclusive_end, more_stats;
        IBD_Pool pool;

    public:
        IBD_Collection(double threshold, bool exclusive_end=true,
                bool more_stats=false) : threshold(threshold),
                exclusive_end(exclusive_end), more_stats(more_stats) {};
        ~IBD_Collection();
        void initialize(int num_samples, Genotype_Reader &reader);
        void initializeWithSites(int num_samples, Genotype_Reader &reader);
        void update(Genotype_Reader &reader, std::ostream &output);
        void purge(std::ostream &output);
};
