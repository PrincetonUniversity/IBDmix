#pragma once

#include <iostream>
#include <vector>

#include "IBDmix/IBD_Stack.h"

class Recorder {
 public:
  virtual void writeHeader(std::ostream &output) const = 0;
  virtual void initializeSegment() = 0;
  virtual void record(IBD_Node *node) = 0;
  virtual void report(std::ostream &output) const = 0;
};

class CountRecorder : public Recorder {
 public:
  void writeHeader(std::ostream &output) const override;
  void initializeSegment() override;
  void record(IBD_Node *node) override;
  void report(std::ostream &output) const override;

 private:
  int in_mask;
  int maf_low;
  int maf_high;
  int rec_2_0;
  int rec_0_2;
  int sites;
  int both;
  int positive_lod;
  int negative_lod;
};

class SiteRecorder : public Recorder {
 public:
  void writeHeader(std::ostream &output) const override;
  void initializeSegment() override;
  void record(IBD_Node *node) override;
  void report(std::ostream &output) const override;

 private:
  std::vector<uint64_t> positions;
};

class LODRecorder : public Recorder {
 public:
  void writeHeader(std::ostream &output) const override;
  void initializeSegment() override;
  void record(IBD_Node *node) override;
  void report(std::ostream &output) const override;

 private:
  std::vector<double> LODs;
};
