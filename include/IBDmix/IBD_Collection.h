#pragma once

#include <iostream>
#include <vector>

#include "IBDmix/Genotype_Reader.h"
#include "IBDmix/IBD_Segment.h"

class IBD_Collection {
 public:
  explicit IBD_Collection(double threshold, bool exclusive_end = true)
      : threshold(threshold), exclusive_end(exclusive_end) {}
  void initialize(const Genotype_Reader &reader);
  void update(const Genotype_Reader &reader, std::ostream &output);
  void purge(std::ostream &output);

  enum Recorder { counts, sites, lods };
  void add_recorder(IBD_Collection::Recorder type);
  void writeHeader(std::ostream &strm) const;

 private:
  std::vector<IBD_Segment> IBDs;
  double threshold;
  bool exclusive_end;
  IBD_Pool pool;
};
