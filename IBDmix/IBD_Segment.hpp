#pragma once
#include "IBD_Stack.hpp"

class IBD_Segment{
    private:
        struct IBD_Node *start, *end, *top;
        int chrom;
        double thresh;
        char *output;
        void add_node(struct IBD_Node *new_node);

    public:
        IBD_Segment(const char *name, int chromosome, double threshold);
        ~IBD_Segment();
        void add_lod(unsigned long int position, double lod, char * output);
        void display(void);
        char *name;
        int length(void);
};
