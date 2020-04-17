#include "IBDmix/IBD_Collection.h"

IBD_Collection::~IBD_Collection(){
    for(auto &ibd : IBDs)
        ibd.reset();
}

void IBD_Collection::initialize(int num_samples, Genotype_Reader &reader){
    IBDs.reserve(num_samples);
    for(auto & sample : reader.get_samples())
        IBDs.push_back(
            std::unique_ptr<IBD_Segment>(
                new IBD_Segment(sample, threshold, &pool,
                    exclusive_end, more_stats)));
}

void IBD_Collection::initializeWithSites(int num_samples, Genotype_Reader &reader){
    IBDs.reserve(num_samples);
    for(auto & sample : reader.get_samples())
        IBDs.push_back(
            std::unique_ptr<IBD_Segment>(
                new IBD_Segment_Sites(sample, threshold, &pool, exclusive_end, more_stats)));
}

void IBD_Collection::update(Genotype_Reader &reader, std::ostream &output){
    for(int i = 0; i < IBDs.size(); i++){
        IBDs[i]->add_lod(reader.chromosome,
                reader.position,
                reader.lod_scores[i],
                reader.line_filtering | reader.recover_type[i],
                output);
    }
}

void IBD_Collection::purge(std::ostream &output){
    for(int i = 0; i < IBDs.size(); i++)
        IBDs[i]->purge(output);
}
