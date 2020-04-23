#include "IBDmix/Segment_Recorders.h"

void CountRecorder::writeHeader(std::ostream &output) const{
    output <<
        "\tsites\tpositive_lods\tmask_and_maf\tin_mask\t"
        "maf_low\tmaf_high\trec_2_0\trec_0_2";
}

void CountRecorder::initializeSegment(){
    positive_lod = both = sites = in_mask =
        maf_low = maf_high = rec_2_0 = rec_0_2 = 0;
}

void CountRecorder::record(IBD_Node *node){
    unsigned char bitmask = node->bitmask;
    if((bitmask & IN_MASK) && ((bitmask & MAF_LOW) || (bitmask & MAF_HIGH))) ++both;
    if((bitmask & IN_MASK) && !(bitmask & MAF_LOW) && !(bitmask &MAF_HIGH)) ++in_mask;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_LOW)) ++maf_low;
    if(!(bitmask & IN_MASK) && (bitmask & MAF_HIGH)) ++maf_high;
    if(bitmask & RECOVER_2_0) ++rec_2_0;
    if(bitmask & RECOVER_0_2) ++rec_0_2;
    if(node->lod > 0) ++positive_lod;
    ++sites;
}

void CountRecorder::report(std::ostream &output) const{
        output << '\t' << sites << '\t'
            << positive_lod << '\t'
            << both << '\t'
            << in_mask << '\t'
            << maf_low << '\t'
            << maf_high << '\t'
            << rec_2_0 << '\t'
            << rec_0_2;
}

void SiteRecorder::writeHeader(std::ostream &output) const{
    output << "\tSNPs";
}

void SiteRecorder::initializeSegment(){
    positions.clear();
}

void SiteRecorder::record(IBD_Node *node){
    if(node->lod > 0)
        positions.push_back(node->position);
}

void SiteRecorder::report(std::ostream &output) const{
    output << '\t';
    for(auto pos = positions.begin(); pos != positions.end(); ++pos){
        if(pos != positions.begin())
            output << ',';
        output << *pos;
    }
}

void LODRecorder::writeHeader(std::ostream &output) const{
    output << "\tLODs";
}

void LODRecorder::initializeSegment(){
    LODs.clear();
}

void LODRecorder::record(IBD_Node *node){
    if(node->lod > 0)
        LODs.push_back(node->lod);
}

void LODRecorder::report(std::ostream &output) const{
    int prec = output.precision();
    output << '\t' << std::setprecision(4);
    for(auto lod = LODs.begin(); lod != LODs.end(); ++lod){
        if(lod != LODs.begin())
            output << ',';
        output << *lod;
    }
    output.precision(prec);  // reset formatting
}
