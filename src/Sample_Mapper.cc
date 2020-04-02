#include "IBDmix/Sample_Mapper.h"

void Sample_Mapper::find_archaic(std::string &archaic){
    if(archaic == ""){
        archaic_index = 0;
        return;
    }

    auto it = std::find(samples.begin(), samples.end(), archaic);
    if(it == samples.end()){
        throw std::invalid_argument("Unable to find archaic '"
                + archaic + '\'');
    }
    archaic_index = std::distance(samples.begin(), it);
}

void Sample_Mapper::map(std::vector<std::string> requested_samples){
    // set the mapping from sample to its index in the genotype file line
    if(requested_samples.empty()){
        //remove archaic index
        samples.erase(samples.begin()+archaic_index);
        sample_to_index.reserve(samples.size());
        // set map to range, removing the archaic index
        int index = 0;
        for(int i = 0; i < samples.size(); ++i, ++index){
            if(index == archaic_index)
                ++index;
            sample_to_index.push_back(index);
        }
    }
    else{
        // ignore archaic index, map to match in buffer
        sample_to_index.reserve(requested_samples.size());
        for(std::string &sample : requested_samples){
            auto it = std::find(samples.begin(), samples.end(), sample);
            if(it == samples.end()){
                throw std::invalid_argument("Unable to find sample '"
                        + sample + '\'');
            }
            sample_to_index.push_back(std::distance(samples.begin(), it));
        }
        // set samples to requested samples, will remove archaic
        samples = requested_samples;
    }
}

int Sample_Mapper::initialize(std::istream &genotype,
        std::istream &requested_samples, std::string archaic){
    // clear any previous results
    sample_to_index.clear();
    samples.clear();

    // read first line of genotype file to populate samples
    std::string line;
    std::string token;
    if(std::getline(genotype, line)){
        std::istringstream iss(line);
        iss >> token;  // chrom
        iss >> token;  // pos
        iss >> token;  // ref
        iss >> token;  // alt
        while(iss >> token)
            samples.emplace_back(token);
    }
    if(samples.empty())
        throw std::invalid_argument("Unable to find samples. Invalid genotype file");

    // if samples is provided, populate requested samples
    std::vector<std::string> requested;
    while(requested_samples >> token)
        requested.emplace_back(token);

    find_archaic(archaic);
    map(requested);
    return samples.size();
}
