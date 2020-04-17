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

    protected:
        virtual void initialize_stats();
        virtual void update_stats(IBD_Node *node);
        virtual void report_stats(std::ostream &output);

    public:
        IBD_Segment(std::string name, double threshold, IBD_Pool *pool,
                bool exclusive_end=true, bool more_stats=false);
        ~IBD_Segment();
        void add_lod(int chromosome, unsigned long int position,
                double lod, unsigned char bitmask, std::ostream &output);
        void purge(std::ostream &output);
        int size();
        void write(std::ostream &strm) const;
};

class IBD_Segment_Sites : public IBD_Segment{
    using IBD_Segment::IBD_Segment;

    private:
        std::vector<unsigned long int> positions;

    protected:
        void initialize_stats();
        void update_stats(IBD_Node *node);
        void report_stats(std::ostream &output);
};


std::ostream& operator<<(std::ostream &strm, const IBD_Segment &segment);
