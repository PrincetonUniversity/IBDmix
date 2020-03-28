#pragma once
#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Genotype_Reader.h"
#include <iostream>
#include <string.h>
#include <limits>

class IBD_Segment{
    private:
        struct IBD_Node *start, *end, *top;
        char *name;
        int chrom;
        double thresh;
        void add_node(struct IBD_Node *new_node, std::ostream &output);
        int in_mask, maf_low, maf_high, rec_2_0, rec_0_2, sites, both;
        void update_counts(unsigned char bitmask);
        bool exclusive_end;
        bool more_stats;

    public:
        IBD_Segment(const char *name, double threshold,
                bool exclusive_end=true, bool more_stats=false);
        ~IBD_Segment();
        void add_lod(int chromosome, unsigned long int position,
                double lod, std::ostream &output, unsigned char bitmask);
        void purge(std::ostream &output);
        int length(void);
        void display(void);
};
