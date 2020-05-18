#include "IBDmix/IBD_Segment.h"

#include <limits>

#include "IBDmix/Genotype_Reader.h"

IBD_Segment::IBD_Segment(std::string segment_name, double threshold,
                         IBD_Pool *pool, bool exclusive_end)
    : name(segment_name),
      threshold(threshold),
      pool(pool),
      exclusive_end(exclusive_end) {}

IBD_Segment::~IBD_Segment() {
  pool->reclaim_stack(&segment);
  for (auto &recorder : recorders) recorder.reset();
}

void IBD_Segment::add_recorder(std::shared_ptr<Recorder> recorder) {
  recorders.push_back(recorder);
}

void IBD_Segment::add_lod(std::string chromosome, uint64_t position, double lod,
                          unsigned char bitmask, std::ostream &output) {
  if (chromosome != "") this->chromosome = chromosome;
  // ignore negative lod as first entry
  if (segment.empty() && lod < 0) {
    return;
  }

  add_node(pool->get_node(position, lod, bitmask), output);
}

void IBD_Segment::purge(std::ostream &output) {
  // -1 and 0s are placeholders, the -inf forces segment to pop all
  add_lod("", 0, -std::numeric_limits<double>::infinity(), 0, output);
}

void IBD_Segment::add_node(IBD_Node *node, std::ostream &output) {
  if (segment.empty() && node->lod < 0) {
    pool->reclaim_node(node);
    return;
  }
  segment.push(node);

  // first entry, reset counts
  if (segment.isSingleton()) {
    initialize_stats();
    update_stats(node);
    return;
  }

  if (segment.topIsNewMax()) {
    // add all nodes from top to end
    if (!recorders.empty()) update_stats_recursive(segment.getTop());
    segment.setEnd();
    pool->reclaim_segment(&segment);
  }

  if (segment.reachedMax()) {
    // write output
    if (segment.endLod() >= threshold) {
      uint64_t pos = segment.endPosition();
      if (exclusive_end && !segment.isEnd(segment.getTop())) {
        // find previous node 'above' end
        const IBD_Node *ptr = segment.getTop();
        while (!segment.isEnd(ptr->next)) ptr = ptr->next;

        if (ptr->lod != -std::numeric_limits<double>::infinity())
          pos = ptr->position;
      }
      output << name << '\t' << chromosome << '\t' << segment.startPosition()
             << '\t' << pos << '\t' << segment.endLod();
      report_stats(output);
      output << '\n';
    }
    IBD_Stack unprocessed = segment.getUnprocessed();

    pool->reclaim_stack(&segment);

    while (!unprocessed.empty()) {
      add_node(unprocessed.pop(), output);
    }
  }
}

void IBD_Segment::initialize_stats() {
  for (auto &recorder : recorders) recorder->initializeSegment();
}

void IBD_Segment::update_stats_recursive(const IBD_Node *node) {
  // since the list is linked in decreasing order, need to traverse in
  // reverse via recursion
  if (node == nullptr || segment.isEnd(node)) return;
  update_stats_recursive(node->next);
  update_stats(node);
}

void IBD_Segment::update_stats(const IBD_Node *node) {
  for (auto &recorder : recorders) recorder->record(node);
}

void IBD_Segment::report_stats(std::ostream &output) {
  for (auto &recorder : recorders) recorder->report(output);
}

void IBD_Segment::writeHeader(std::ostream &strm) const {
  for (auto &recorder : recorders) recorder->writeHeader(strm);
}

void IBD_Segment::write(std::ostream &strm) const {
  strm << "--- " << name << " ---\n";
  segment.write(strm);
}

int IBD_Segment::size(void) const { return segment.size(); }

std::ostream &operator<<(std::ostream &strm, const IBD_Segment &segment) {
  segment.write(strm);
  return strm;
}
