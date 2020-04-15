#pragma once
#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Genotype_Reader.h"
#include <iostream>
#include <string>
#include <limits>

class IBD_Segment{
    private:
        std::string name;
        IBD_Stack segment;
        IBD_Pool * pool;
        int chrom;
        double thresh;
        int in_mask, maf_low, maf_high, rec_2_0, rec_0_2, sites, both;
        bool exclusive_end;
        bool more_stats;

        void add_node(IBD_Node *node, std::ostream &output);
        void update_counts(unsigned char bitmask);

    public:
        IBD_Segment(std::string name, double threshold, IBD_Pool *pool,
                bool exclusive_end=true, bool more_stats=false);
        ~IBD_Segment();
        void add_lod(int chromosome, unsigned long int position,
                double lod, unsigned char bitmask, std::ostream &output);
        void purge(std::ostream &output);
        int size();
        void display();
};
