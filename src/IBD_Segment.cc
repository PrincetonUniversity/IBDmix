#include "IBDmix/IBD_Segment.h"

#include <limits>

#include "IBDmix/Genotype_Reader.h"

IBD_Segment::IBD_Segment(std::string segment_name, double threshold,
        IBD_Pool *pool, bool exclusive_end) :
    name(segment_name), thresh(threshold), pool(pool),
    exclusive_end(exclusive_end) {}

IBD_Segment::~IBD_Segment() {
    pool->reclaim_all(&segment.top);
    for (auto & recorder : recorders)
        recorder.reset();
}

void IBD_Segment::add_recorder(std::shared_ptr<Recorder> recorder) {
    recorders.push_back(recorder);
}


void IBD_Segment::add_lod(int chromosome, uint64_t position,
        double lod, unsigned char bitmask, std::ostream &output) {
    if (chromosome > 0) chrom = chromosome;
    // ignore negative lod as first entry
    if (segment.empty() && lod < 0) {
        return;
    }

    add_node(pool->get_node(position, lod, bitmask), output);
}

void IBD_Segment::purge(std::ostream &output) {
    // -1 and 0s are placeholders, the -inf forces segment to pop all
    add_lod(-1, 0, -std::numeric_limits<double>::infinity(), 0, output);
}

void IBD_Segment::add_node(IBD_Node *node, std::ostream &output) {
    if (segment.empty() && node->lod < 0) {
        pool->reclaim_node(node);
        return;
    }
    segment.push(node);

    // first entry, reset counts
    if (segment.top->next == nullptr) {
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
    if (segment.top->cumulative_lod >= segment.end->cumulative_lod) {
        // add all nodes from top to end
        if (!recorders.empty())
            update_stats_recursive(segment.top);
        segment.end = segment.top;
        pool->reclaim_between(segment.end, segment.start);
    }

    if (segment.top->cumulative_lod < 0) {
        // write output
        if (segment.end->cumulative_lod >= thresh) {
            uint64_t pos = segment.end->position;
            if (exclusive_end && segment.end != segment.top) {
                // find previous node 'above' top
                struct IBD_Node * ptr = segment.top;
                for (; ptr->next != segment.end; ptr=ptr->next) {}
                if (ptr->lod != -std::numeric_limits<double>::infinity())
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
        IBD_Stack reversed(segment.end->next);
        segment.end->next = nullptr;
        pool->reclaim_all(&segment.top);
        // reset member variables to process reversed
        segment.top = segment.start = segment.end = nullptr;
        while (!reversed.empty()) {
            add_node(reversed.pop(), output);
        }
    }
}

void IBD_Segment::initialize_stats() {
    for ( auto &recorder : recorders)
        recorder->initializeSegment();
}

void IBD_Segment::update_stats_recursive(IBD_Node *node) {
    // since the list is linked in decreasing order, need to traverse in
    // reverse via recursion
    if (node == segment.end || node == nullptr)
        return;
    update_stats_recursive(node->next);
    update_stats(node);
}

void IBD_Segment::update_stats(IBD_Node *node) {
    for ( auto &recorder : recorders)
        recorder->record(node);
}

void IBD_Segment::report_stats(std::ostream &output) {
    for ( auto &recorder : recorders)
        recorder->report(output);
}

void IBD_Segment::writeHeader(std::ostream &strm) const {
    for ( auto &recorder : recorders)
        recorder->writeHeader(strm);
}

void IBD_Segment::write(std::ostream &strm) const {
    strm << "--- " << name << " ---\n";
    for (struct IBD_Node* ptr = segment.top; ptr != nullptr; ptr = ptr->next) {
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

int IBD_Segment::size(void) {
    return segment.size();
}

std::ostream& operator<<(std::ostream &strm, const IBD_Segment &segment) {
    segment.write(strm);
    return strm;
}
