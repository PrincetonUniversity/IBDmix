#include "IBDmix/IBD_Segment.h"

IBD_Segment::IBD_Segment(std::string segment_name, double threshold,
        IBD_Pool *pool, bool exclusive_end, bool more_stats) :
    name(segment_name), thresh(threshold), pool(pool),
    exclusive_end(exclusive_end), more_stats(more_stats)
{ }

IBD_Segment::~IBD_Segment(){
    pool->reclaim_all(segment.top);
}

void IBD_Segment::add_lod(int chromosome, unsigned long int position,
        double lod, unsigned char bitmask, std::ostream &output){
    if(chromosome > 0) chrom = chromosome;
    // ignore negative lod as first entry
    if(segment.empty() && lod < 0){
        return;
    }

    add_node(pool->get_node(position, lod, bitmask), output);
}

void IBD_Segment::purge(std::ostream &output){
    // -1 and 0s are placeholders, the -inf forces segment to pop all
    add_lod(-1, 0, -std::numeric_limits<double>::infinity(), 0, output);
}

void IBD_Segment::add_node(IBD_Node *node, std::ostream &output){
    if(segment.empty() && node->lod < 0){
        pool->reclaim_node(node);
        return;
    }
    segment.push(node);

    // first entry, reset counts
    if(segment.top->next == nullptr){
        both = sites = in_mask = maf_low = maf_high = rec_2_0 = rec_0_2 = 0;
        update_counts(node->bitmask);
        segment.start = segment.end = segment.top;
        node->cumulative_lod = node->lod;
        return;
    }
    // cumulative lod defaults to lod
    segment.top->cumulative_lod = segment.top->next->cumulative_lod +
        segment.top->lod;

    // new max, collapse to start
    if(segment.top->cumulative_lod >= segment.end->cumulative_lod){
        // TODO does this need to be all between?
        update_counts(node->bitmask);
        segment.end = segment.top;
        pool->reclaim_between(segment.end, segment.start);
    }

    if(segment.top->cumulative_lod < 0){
        // write output
        if(segment.end->cumulative_lod >= thresh){
            unsigned long pos = segment.end->position;
            if(exclusive_end && segment.end != segment.top){
                //find previous node 'above' top
                struct IBD_Node * ptr = segment.top;
                for(; ptr->next != segment.end; ptr=ptr->next);
                if(ptr->lod != -std::numeric_limits<double>::infinity())
                    pos = ptr->position;
            }
            output << name << '\t'
                << chrom << '\t'
                << segment.start->position << '\t'
                << pos << '\t'
                << segment.end->cumulative_lod;
            if (more_stats == true)
                output << '\t' << sites << '\t'
                    << both << '\t'
                    << in_mask << '\t'
                    << maf_low << '\t'
                    << maf_high << '\t'
                    << rec_2_0 << '\t'
                    << rec_0_2;
            output << '\n';
        }
        // reverse list to reprocess remaining nodes
        segment.reverse();
        // skip from top to end
        IBD_Stack reversed = segment.end->next;
        segment.end->next = nullptr;
        pool->reclaim_all(segment.top);
        // reset member variables to process reversed
        segment.top = segment.start = segment.end = nullptr;
        while(!reversed.empty()){
            add_node(reversed.pop(), output);
        }
    }
}

void IBD_Segment::update_counts(unsigned char bitmask){
    if (more_stats == false)
        return;
    if((bitmask & IN_MASK) && ((bitmask & MAF_LOW) || (bitmask & MAF_HIGH))) both++;
    if((bitmask & IN_MASK) && !(bitmask & MAF_LOW) && !(bitmask &MAF_HIGH)) in_mask++;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_LOW)) maf_low++;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_HIGH)) maf_high++;
    if(bitmask & RECOVER_2_0) rec_2_0++;
    if(bitmask & RECOVER_0_2) rec_0_2++;
    sites++;
}

int IBD_Segment::size(void){
    return segment.size();
}

void IBD_Segment::write(std::ostream &strm) const{
    strm << "--- " << name << " ---\n";
    for(struct IBD_Node* ptr = segment.top; ptr != nullptr; ptr = ptr->next){
        strm << ptr->position << "\t"
            << ptr->lod << "\t"
            << ptr->cumulative_lod;
        if (ptr == segment.top)
            strm << " <- top";
        if (ptr == segment.start)
            strm << " <- start";
        if (ptr == segment.end)
            strm << " <- end";
        strm << "\n";
    }
}

std::ostream& operator<<(std::ostream &strm, const IBD_Segment &segment){
    segment.write(strm);
    return strm;
}
