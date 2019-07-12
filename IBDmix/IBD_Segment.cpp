#include <iostream>
#include <string.h>
#include <limits>
#include "IBD_Segment.hpp"
#include "IBD_Stack.hpp"

IBD_Segment::IBD_Segment(const char *segment_name, double threshold){
    // NOTE!! start and end refer to the segment, not the stack
    // The stack grows nearest end, start is nearest to bottom
    start = end = top = nullptr;
    name = new char[strlen(segment_name)+1];
    strcpy(name, segment_name);
    thresh = threshold;
}

IBD_Segment::~IBD_Segment(){
    delete[] name;
}

void IBD_Segment::add_lod(int chromosome, unsigned long int position,
        double lod, FILE *output){
    chrom = chromosome;
    add_node(get_node(position, lod), output);
}

void IBD_Segment::purge(FILE *output){
    add_node(get_node(-1, -std::numeric_limits<double>::infinity()),
            output);
}

void IBD_Segment::add_node(struct IBD_Node *new_node, FILE * output){
    // ignore negative lod as first entry
    if(top == nullptr && new_node->lod < 0){
        reclaim_node(new_node);
        return;
    }

    push(top, new_node);

    // first entry, nothing else to do
    if(top->next == nullptr){
        start = end = top;
        top->cumulative_lod = top->lod;
        return;
    }
    top->cumulative_lod = top->next->cumulative_lod + top->lod;

    // new max, collapse to start
    if(top->cumulative_lod >= end->cumulative_lod){
        end = top;
        reclaim_between(end, start);
    }

    if(top->cumulative_lod < 0){
        // write output
        if(end->cumulative_lod >= thresh){
            // TODO this is slow and bad to match legacy, probably change
            unsigned long int pos = end->position;
            if(end != top){
                //find previous node 'above' top
                struct IBD_Node * ptr = top;
                for(; ptr->next != end; ptr=ptr->next);
                if(ptr->lod != -std::numeric_limits<double>::infinity())
                    pos = ptr->position;
            }
            fprintf(output, "%s\t%d\t%lu\t%lu\t%g\n",
                    name, chrom, start->position,
                    pos, end->cumulative_lod);
        }
        // reverse list to reprocess remaining nodes
        top = reverse(top);
        // skip from top to end
        struct IBD_Node *reversed = end->next;
        end->next = nullptr;
        reclaim_all(top);
        // reset member variables to process reversed
        top = start = end = nullptr;
        while(reversed != nullptr)
            add_node(pop(reversed), output);
    }
}

int IBD_Segment::length(void){
    return stack_length(top);
}

void IBD_Segment::display(void){
    std::cout << "--- " << name << " ---\n";
    for(struct IBD_Node* ptr = top; ptr != nullptr; ptr = ptr->next){
        std::cout << ptr->position << "\t"
            << ptr-> lod << "\t"
            << ptr-> cumulative_lod;
        if (ptr == top)
            std::cout << " <- top";
        if (ptr == start)
            std::cout << " <- start";
        if (ptr == end)
            std::cout << " <- end";
        std::cout << "\n";
    }
}
