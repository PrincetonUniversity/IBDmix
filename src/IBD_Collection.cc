#include <stdlib.h>
#include <string.h>
#include "IBDmix/IBD_Collection.h"
#include "IBDmix/Genotype_Reader.h"

IBD_Collection::IBD_Collection(){
    num_samples = 0;
}

void IBD_Collection::initialize(int num_samples, double threshold,
        Genotype_Reader &reader, bool exclusive_end, bool more_stats){
    this->num_samples = num_samples;
    IBDs.reserve(num_samples);
    char *sample;
    for(int i = 0; i < num_samples; i++){
        reader.yield_sample(sample, i);
        IBDs.emplace_back(sample, threshold, exclusive_end, more_stats);
    }
}

void IBD_Collection::update(Genotype_Reader &reader, FILE * output){
    for(int i = 0; i < num_samples; i++){
        IBDs[i].add_lod(reader.chromosome,
                reader.position,
                reader.lod_scores[i],
                output,
                reader.line_filtering | reader.recover_type[i]);
    }
}

void IBD_Collection::purge(FILE * output){
    for(int i = 0; i < num_samples; i++)
        IBDs[i].purge(output);
}
