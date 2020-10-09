#include "IBDmix/IBD_Collection.h"

void IBD_Collection::initialize(const Genotype_Reader &reader) {
  IBDs.reserve(reader.get_samples().size());
  for (auto &sample : reader.get_samples())
    IBDs.emplace_back(sample, threshold, &pool, exclusive_end);
}

void IBD_Collection::update(const Genotype_Reader &reader,
                            std::ostream &output) {
  for (unsigned int i = 0; i < IBDs.size(); i++) {
    IBDs[i].add_lod(reader.getChromosome(), reader.getPosition(),
                    reader.getLodScore(i),
                    reader.getLineFilter() | reader.getRecoverType(i), output);
  }
}

void IBD_Collection::purge(std::ostream &output) {
  for (unsigned int i = 0; i < IBDs.size(); i++) IBDs[i].purge(output);
}

void IBD_Collection::add_recorder(IBD_Collection::Recorder type) {
  switch (type) {
    case IBD_Collection::Recorder::counts:
      for (auto &ibd : IBDs)
        ibd.add_recorder(std::make_shared<CountRecorder>());
      break;
    case IBD_Collection::Recorder::sites:
      for (auto &ibd : IBDs) ibd.add_recorder(std::make_shared<SiteRecorder>());
      break;
    case IBD_Collection::Recorder::lods:
      for (auto &ibd : IBDs) ibd.add_recorder(std::make_shared<LODRecorder>());
      break;
  }
}

void IBD_Collection::writeHeader(std::ostream &strm) const {
  IBDs[0].writeHeader(strm);
}
