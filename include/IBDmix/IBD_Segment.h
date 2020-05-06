#ifndef INCLUDE_IBDMIX_IBD_SEGMENT_H_
#define INCLUDE_IBDMIX_IBD_SEGMENT_H_

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Segment_Recorders.h"

class IBD_Segment{
 private:
    std::string name;
    double thresh;
    IBD_Stack segment;
    IBD_Pool * pool;
    std::vector<std::shared_ptr<Recorder>> recorders;
    int chrom;
    bool exclusive_end;

    void add_node(IBD_Node *node, std::ostream &output);
    void initialize_stats();
    void update_stats_recursive(IBD_Node *node);
    void update_stats(IBD_Node *node);
    void report_stats(std::ostream &output);

 public:
    IBD_Segment(std::string name, double threshold, IBD_Pool *pool,
            bool exclusive_end = true);
    ~IBD_Segment();
    void add_lod(int chromosome, uint64_t position,
            double lod, unsigned char bitmask, std::ostream &output);
    void add_recorder(std::shared_ptr<Recorder> recorder);
    void purge(std::ostream &output);
    int size();
    void write(std::ostream &strm) const;
    void writeHeader(std::ostream &strm) const;
};

std::ostream& operator<<(std::ostream &strm, const IBD_Segment &segment);

#endif  // INCLUDE_IBDMIX_IBD_SEGMENT_H_
