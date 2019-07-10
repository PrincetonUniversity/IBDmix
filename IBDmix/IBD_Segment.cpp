#include <iostream>
#include <string.h>
#include <stdio.h>
#include "IBD_Segment.hpp"
#include "IBD_Stack.hpp"

IBD_Segment::IBD_Segment(const char *segment_name, int chromosome, double threshold){
    // NOTE!! start and end refer to the segment, not the stack
    // The stack grows nearest end, start is nearest to bottom
    start = end = top = nullptr;
    name = new char[strlen(segment_name)];
    strcpy(name, segment_name);
    chrom = chromosome;
    thresh = threshold;
}

IBD_Segment::~IBD_Segment(){
    delete name;
    reclaim_all(top);
}

void IBD_Segment::add_lod(unsigned long int position, double lod, char *output){
    // note, output needs to be a member function to keep track of where
    // the current output is between calls of add_node
    this->output = output;
    add_node(get_node(position, lod));
}

void IBD_Segment::add_node(struct IBD_Node *new_node){
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
        // write output, move pointer by size
        if(end->cumulative_lod >= thresh)
            output += sprintf(output, "%s\t%d\t%lu\t%lu\t%g\n",
                    name, chrom, start->position,
                    end->position, end->cumulative_lod);
        // reverse list to reprocess remaining nodes
        top = reverse(top);
        // skip from top to end
        struct IBD_Node *reversed = end->next;
        end->next = nullptr;
        reclaim_all(top);
        // reset member variables to process reversed
        top = start = end = nullptr;
        while(reversed != nullptr)
            add_node(pop(reversed));
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
