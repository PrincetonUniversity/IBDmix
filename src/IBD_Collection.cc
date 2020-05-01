#include "IBDmix/IBD_Collection.h"

void IBD_Collection::initialize(Genotype_Reader &reader){
    IBDs.reserve(reader.get_samples().size());
    for(auto & sample : reader.get_samples())
        IBDs.emplace_back(sample, threshold, &pool, exclusive_end);
}

void IBD_Collection::update(Genotype_Reader &reader, std::ostream &output){
    for(int i = 0; i < IBDs.size(); i++){
        IBDs[i].add_lod(reader.chromosome,
                reader.position,
                reader.lod_scores[i],
                reader.line_filtering | reader.recover_type[i],
                output);
    }
}

void IBD_Collection::purge(std::ostream &output){
    for(int i = 0; i < IBDs.size(); i++)
        IBDs[i].purge(output);
}

void IBD_Collection::add_recorder(IBD_Collection::Recorder type){
    switch(type){
        case IBD_Collection::Recorder::counts:
            for(auto &ibd : IBDs)
                ibd.add_recorder(std::make_shared<CountRecorder>());
            break;
        case IBD_Collection::Recorder::sites:
            for(auto &ibd : IBDs)
                ibd.add_recorder(std::make_shared<SiteRecorder>());
            break;
    }
}

void IBD_Collection::writeHeader(std::ostream &strm) const{
    IBDs[0].writeHeader(strm);
}
