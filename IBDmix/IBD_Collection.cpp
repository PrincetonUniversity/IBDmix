#include <stdlib.h>
#include <string.h>
#include "IBD_Collection.hpp"
#include "Genotype_Reader.hpp"

IBD_Collection::IBD_Collection(){
    num_samples = 0;
}

void IBD_Collection::initialize(int num_samples, double threshold,
        Genotype_Reader &reader){
    this->num_samples = num_samples;
    IBDs.reserve(num_samples);
    char *sample;
    for(int i = 0; i < num_samples; i++){
        reader.yield_sample(sample, i);
        IBDs.emplace_back(sample, threshold);
    }
}

void IBD_Collection::update(Genotype_Reader &reader, FILE * output, FILE * sample){
    fprintf(sample, "%lu\t%f\t%d\t%d\n", reader.position, reader.lod_scores[0],
            reader.line_filtering & IN_MASK ? 1 : 0,
            reader.line_filtering & MAF_LOW ? 1 : 0);
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
