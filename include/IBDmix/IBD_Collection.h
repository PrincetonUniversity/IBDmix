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
        std::vector<IBD_Segment> IBDs;
        double threshold;
        bool exclusive_end;
        IBD_Pool pool;

    public:
        IBD_Collection(double threshold, bool exclusive_end=true) :
            threshold(threshold), exclusive_end(exclusive_end){};
        void initialize(Genotype_Reader &reader);
        void update(Genotype_Reader &reader, std::ostream &output);
        void purge(std::ostream &output);

        enum Recorder { counts, sites, lods };
        void add_recorder(IBD_Collection::Recorder type);
        void writeHeader(std::ostream &strm) const;
};
