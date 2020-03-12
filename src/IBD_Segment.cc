#include <iostream>
#include <string.h>
#include <limits>
#include "IBDmix/IBD_Segment.h"
#include "IBDmix/IBD_Stack.h"

IBD_Segment::IBD_Segment(const char *segment_name, double threshold,
        bool exclusive_end, bool more_stats){
    // NOTE!! start and end refer to the segment, not the stack
    // The stack grows nearest end, start is nearest to bottom
    start = end = top = nullptr;
    name = new char[strlen(segment_name)+1];
    strcpy(name, segment_name);
    thresh = threshold;
    this->exclusive_end = exclusive_end;  // true to match legacy of 'next' position
    this->more_stats = more_stats;
}

IBD_Segment::~IBD_Segment(){
    delete[] name;
}

void IBD_Segment::add_lod(int chromosome, unsigned long int position,
        double lod, FILE *output, unsigned char bitmask){
    chrom = chromosome;
    add_node(get_node(position, lod, bitmask), output);
}

void IBD_Segment::purge(FILE *output){
    add_node(get_node(-1, -std::numeric_limits<double>::infinity()), output);
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
        both = sites = in_mask = maf_low = maf_high = rec_2_0 = rec_0_2 = 0;
        update_counts(new_node->bitmask);
        start = end = top;
        top->cumulative_lod = top->lod;
        return;
    }
    top->cumulative_lod = top->next->cumulative_lod + top->lod;

    // new max, collapse to start
    if(top->cumulative_lod >= end->cumulative_lod){
        update_counts(new_node->bitmask);
        end = top;
        reclaim_between(end, start);
    }

    if(top->cumulative_lod < 0){
        // write output
        if(end->cumulative_lod >= thresh){
            unsigned long int pos = end->position;
            if(exclusive_end && end != top){
                //find previous node 'above' top
                struct IBD_Node * ptr = top;
                for(; ptr->next != end; ptr=ptr->next);
                if(ptr->lod != -std::numeric_limits<double>::infinity())
                    pos = ptr->position;
            }
            fprintf(output, "%s\t%d\t%lu\t%lu\t%g",
                    name, chrom, start->position,
                    pos, end->cumulative_lod);
            if (more_stats == true)
                fprintf(output, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                        sites, both, in_mask, maf_low,
                        maf_high, rec_2_0, rec_0_2);
            fprintf(output, "\n");
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

void IBD_Segment::update_counts(unsigned char bitmask){
    if (more_stats == false)
        return;
    if((bitmask & IN_MASK) && ((bitmask & MAF_LOW) || (bitmask & MAF_HIGH))) both++;
    if((bitmask & IN_MASK) && !(bitmask & MAF_LOW)) in_mask++;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_LOW)) maf_low++;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_HIGH)) maf_high++;
    if(bitmask & RECOVER_2_0) rec_2_0++;
    if(bitmask & RECOVER_0_2) rec_0_2++;
    sites++;
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
