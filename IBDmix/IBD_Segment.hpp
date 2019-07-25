#pragma once
#include "IBD_Stack.hpp"
#include <stdio.h>

class IBD_Segment{
    private:
        struct IBD_Node *start, *end, *top;
        char *name;
        int chrom;
        double thresh;
        void add_node(struct IBD_Node *new_node, FILE *output);

    public:
        IBD_Segment(const char *name, double threshold);
        ~IBD_Segment();
        void add_lod(int chromosome, unsigned long int position,
                double lod, FILE * output);
        void purge(FILE * output);
        void display(void);
        int length(void);
};
