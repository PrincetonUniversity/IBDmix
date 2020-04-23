#pragma once
#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Genotype_Reader.h"
#include "IBDmix/Segment_Recorders.h"
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <memory>

class IBD_Segment{
    private:
        std::string name;
        IBD_Stack segment;
        IBD_Pool * pool;
        std::vector<std::shared_ptr<Recorder>> recorders;
        int chrom;
        double thresh;
        bool exclusive_end;

        void add_node(IBD_Node *node, std::ostream &output);
        void initialize_stats();
        void update_stats(IBD_Node *node);
        void report_stats(std::ostream &output);

    public:
        IBD_Segment(std::string name, double threshold, IBD_Pool *pool,
                bool exclusive_end=true);
        ~IBD_Segment();
        void add_lod(int chromosome, unsigned long int position,
                double lod, unsigned char bitmask, std::ostream &output);
        void add_recorder(std::shared_ptr<Recorder> recorder);
        void purge(std::ostream &output);
        int size();
        void write(std::ostream &strm) const;
        void writeHeader(std::ostream &strm) const;
};

std::ostream& operator<<(std::ostream &strm, const IBD_Segment &segment);
