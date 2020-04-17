#include "IBDmix/IBD_Segment.h"

IBD_Segment::IBD_Segment(std::string segment_name, double threshold,
        IBD_Pool *pool, bool exclusive_end, bool more_stats) :
    name(segment_name), thresh(threshold), pool(pool),
    exclusive_end(exclusive_end), more_stats(more_stats){
        both = sites = in_mask = maf_low = maf_high = rec_2_0 = rec_0_2 = 0;
    }

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
        initialize_stats();
        update_stats(node);
        segment.start = segment.end = segment.top;
        node->cumulative_lod = node->lod;
        return;
    }
    // cumulative lod defaults to lod
    segment.top->cumulative_lod = segment.top->next->cumulative_lod +
        segment.top->lod;

    // new max, collapse to start
    if(segment.top->cumulative_lod >= segment.end->cumulative_lod){
        update_stats(node);
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
            report_stats(output);
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

void IBD_Segment::initialize_stats(){
    both = sites = in_mask = maf_low = maf_high = rec_2_0 = rec_0_2 = 0;
}

void IBD_Segment::update_stats(IBD_Node *node){
    if (more_stats == false)
        return;
    unsigned char bitmask = node->bitmask;
    if((bitmask & IN_MASK) && ((bitmask & MAF_LOW) || (bitmask & MAF_HIGH))) both++;
    if((bitmask & IN_MASK) && !(bitmask & MAF_LOW) && !(bitmask &MAF_HIGH)) in_mask++;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_LOW)) maf_low++;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_HIGH)) maf_high++;
    if(bitmask & RECOVER_2_0) rec_2_0++;
    if(bitmask & RECOVER_0_2) rec_0_2++;
    sites++;
}

void IBD_Segment::report_stats(std::ostream &output){
    if (more_stats == true){
        output << '\t' << sites << '\t'
            << both << '\t'
            << in_mask << '\t'
            << maf_low << '\t'
            << maf_high << '\t'
            << rec_2_0 << '\t'
            << rec_0_2;
    }
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

void IBD_Segment_Sites::initialize_stats(){
    IBD_Segment::initialize_stats();
    positions.clear();
}

void IBD_Segment_Sites::update_stats(IBD_Node *node){
    IBD_Segment::update_stats(node);
    if(node->lod > 0)
        positions.push_back(node->position);
}

void IBD_Segment_Sites::report_stats(std::ostream &output){
    IBD_Segment::report_stats(output);
    output << '\t';
    for(const auto& pos : positions){
        if(pos != positions[0])
            output << ',';
        output << pos;
    }
}

int IBD_Segment::size(void){
    return segment.size();
}

std::ostream& operator<<(std::ostream &strm, const IBD_Segment &segment){
    segment.write(strm);
    return strm;
}
