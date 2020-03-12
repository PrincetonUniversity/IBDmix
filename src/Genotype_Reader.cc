#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include "IBDmix/Genotype_Reader.h"

bool is_break(const char * chr);

Genotype_Reader::Genotype_Reader(FILE * genotype, FILE * mask,
        double archaic_error, double modern_error_max,
        double modern_error_proportion, double minesp,
        int minor_allele_cutoff){
    this->genotype = genotype;
    this->mask = mask;
    this->archaic_error = archaic_error;
    this->modern_error_max = modern_error_max;
    this->modern_error_proportion = modern_error_proportion;
    this->minesp = minesp;
    this->minor_allele_cutoff = minor_allele_cutoff;

    buffer = NULL;
    buf_size = 0;
    samples = NULL;
    sample_to_index = nullptr;
    lod_scores = nullptr;
    recover_type = nullptr;
    lod_cache = new double[3];
    mask_chromosome = -1; // force read on first check


    both = num_lines = count_in_mask = fail_maf = count_recovered = 0;
    frac_rec = 0;
}

Genotype_Reader::~Genotype_Reader(){
    if(buffer)
        free(buffer);
    if(samples && samples != buffer)  // can equal if buffer is unused
        free(samples);
    if(sample_to_index != nullptr)
        delete []sample_to_index;
    if(lod_scores != nullptr)
        delete []lod_scores;
    if(recover_type != nullptr)
        delete []recover_type;
    delete []lod_cache;
}

int Genotype_Reader::initialize(FILE * samples, const char * archaic){
    // using samples list and header line, determine number of samples
    // and mapping from position to sample number
    int result = 0;
    //get number of samples
    char * sample_line = NULL;
    size_t sample_len = 0;

    //read genotype header
    fscanf(genotype, "%*s\t%*s\t%*s\t%*s\t");  // chrom, pos, ref, alt
    getline(&buffer, &buf_size, genotype);

    if(samples == nullptr){
        for(char *ptr = buffer; *ptr != '\0'; ptr++)
            if(*ptr == '\t')
                result++;
        result++; // last token has newline instead of tab
        result--; // but remove one for the archaic sample
        this->samples = buffer;
    }
    else{
        size_t read = getdelim(&sample_line, &sample_len, '\0', samples);
        for(char *ptr = sample_line; *ptr != '\0'; ptr++)
            if(*ptr == '\n'){
                *ptr = '\t'; //make tab delimited
                result++;
            }

        // handle case where sample file is not newline terminated
        if (sample_line[read-1] != '\t')
            result++;
        this->samples = sample_line;
    }
    num_samples = result;
    sample_to_index = new int[result];
    find_archaic(archaic, sample_line);
    determine_sample_mapping(sample_line);

    lod_scores = new double[result];
    recover_type = new unsigned char[result];
    return result;
}

void Genotype_Reader::determine_sample_mapping(const char * sample_line){
    // set the mapping from sample to its index in the genotype file line
    if(sample_line == NULL){
        // set map to range, removing the archaic index
        int index = 0;
        for(int i = 0; i < num_samples; i++, index++){
            if(index == archaic_index)
                index++;
            sample_to_index[i] = index;
        }
    }
    else{
        // ignore archaic index, map to match in buffer
        const char * next_sample = sample_line;
        for(int i = 0; i < num_samples; i++){
            if( (sample_to_index[i] = find_token(next_sample, buffer)) == -1){
                std::cout << "Unable to find sample " << i << ": \"";
                for(; *next_sample != '\0' && *next_sample != '\t'; next_sample++)
                    std::cout << next_sample[0];
                std::cout << "\"\n";
                exit(1);
            }
            for(; *next_sample != '\0' && *next_sample != '\t'; next_sample++);
            if(*next_sample == '\t') next_sample++;
        }
    }
}

void Genotype_Reader::find_archaic(const char * archaic,
        const char * sample_line){
    if(archaic == nullptr){
        // set sample map
        archaic_index = 0;
        if (sample_line == nullptr){
            // remove first name from samples
            int len = 0; // +1 for the tab
            char * sp;
            for(sp = samples; *sp != '\t' && *sp != '\0'; sp++)
                len++;
            if(*sp == '\t') len++;
            for(sp = samples; *(sp+len) != '\0'; sp++)
                *sp = *(sp + len);
            *sp = '\0';

        }
        // else keep archaic in samples regardless
    }
    else{
        // find position of archaic
        archaic_index = find_token(archaic, buffer);
        if(archaic_index == -1){
            std::cout << "Unable to find archaic: '" << archaic << "'\n";
            exit(1);
        }
        // remove from samples
        if (sample_line == nullptr){
            char * sp = samples;
            // go to token of archaic
            for(int count = 0; count != archaic_index; count++){
                for(; *sp != '\t'; sp++);
                sp++;
            }
            // copy over string
            int len = strlen(archaic);
            if(*(sp + len) == '\t') len++;  //skip tab if it's next character
            while(*(sp + len) != '\0'){
                *sp = *(sp + len);
                sp++;
            }
            *sp = '\0';
        }
    }
}

bool Genotype_Reader::update(){
    // read next line of input file
    // update the lod_scores array for reading, handling masks
    // return false if the file is read fully
    if(fscanf(genotype, "%i\t%lu\t%*s\t%*s\t", &chromosome, &position) == EOF){
        return false;
    }
    line_filtering = 0;
    getline(&buffer, &buf_size, genotype);
    num_lines++;
    // selected indicates if the line should have its lod calculated
    // set to false if one of the following occurs:
    // - in a masked region
    // - fails to meet allele cutoff
    // If selected is false, lod = 0, unless archaic = (0, 2) and modern = (2, 0)
    bool selected = !in_mask();
    if(!selected){
        count_in_mask++;
        line_filtering |= IN_MASK;
    }
    process_line_buffer(selected);
    return true;
}

bool Genotype_Reader::in_mask(){
    //check if the current chromosome, position is within the mask region
    //assumes both are in sorted order and mask is merged!
    if(mask == nullptr)
        return false;

    while(!feof(mask) && (mask_chromosome < chromosome ||
            (mask_chromosome == chromosome && position > mask_end)))
        if (fscanf(mask, "%i%lu%lu\n", &mask_chromosome, &mask_start, &mask_end) == 0){
            printf("Unable to read mask file; chromosome must be an integer\n");
            exit(1);
        }

    return mask_chromosome == chromosome && position > mask_start &&
        position <= mask_end;

}

void Genotype_Reader::process_line_buffer(bool selected){
    // assume buffer is loaded with tab-separated character in {0, 1, 2, 9}
    // from a genotype file.  Using sample_to_index mapping, fill in
    // the lod_scores array with appropriate values
    // All arrays must be initialized!
    char archaic = buffer[archaic_index*2];  //throughout, *2 to skip tabs
    double allele_frequency = 0;
    bool temp = get_frequency(allele_frequency);
    if(!temp) fail_maf++;
    selected &= temp;
    if(!selected) both++;
    double modern_error = get_modern_error(allele_frequency);

    update_lod_cache(archaic, allele_frequency, modern_error, selected);
    bool recovered = false;
    for(int i = 0; i < num_samples; i++){
        recover_type[i] = 0;
        if(!selected){
            if(archaic == '0' && buffer[sample_to_index[i]*2] == '2'){
                recovered = true;
                frac_rec++;
                recover_type[i] |= RECOVER_0_2;
            }
            if((archaic == '2' && buffer[sample_to_index[i]*2] == '0')){
                recovered = true;
                frac_rec++;
                recover_type[i] |= RECOVER_2_0;
            }
        }
        lod_scores[i] = calculate_lod(buffer[sample_to_index[i]*2]);
    }
    if(recovered) count_recovered++;
}

bool Genotype_Reader::get_frequency(double &frequency){
    // determine the observed frequency of alternative alleles
    // Returns true if enough counts were observed above the cutoff value
    int total_counts=0, alt_counts=0;
    char current;
    bool select = true;
    for(int i = 0; i < num_samples; i++)
        if((current = buffer[sample_to_index[i]*2]) != '9'){
            total_counts += 2;
            alt_counts += current - '0';
        }

    if (alt_counts <= minor_allele_cutoff){ // not enough counts
        select = false;
        line_filtering |= MAF_LOW;
    }
    if (total_counts - alt_counts <= minor_allele_cutoff){ // too many
        select = false;
        line_filtering |= MAF_HIGH;
    }

    if(total_counts == 0)
        frequency = 0;
    else
        frequency = ((double)alt_counts) / total_counts;
    return select;
}

double Genotype_Reader::get_modern_error(double frequency){
    if(frequency > 0.5)  // convert to minor frequency
        frequency = 1 - frequency;
    double prop_error = frequency * modern_error_proportion;

    return prop_error < modern_error_max ? prop_error : modern_error_max;
}

void Genotype_Reader::update_lod_cache(char archaic, double freq_b,
        double modern_error, bool selected){
    // update the lod cache array based on values along genotype file line
    // TODO should minesp be added to more terms? (e.g. 2,2)
    double freq_a = 1 - freq_b;
    double err_0 = 1 - archaic_error * (1 - archaic_error);
    double err_1 = (1 - archaic_error) * (1 - modern_error) + archaic_error * modern_error;
    double err_2 = (1 - archaic_error) * modern_error + archaic_error * (1 - modern_error);

    if(archaic == '0'){
        if(selected){
            lod_cache[0] = log10(err_1 / freq_a / err_0);
            lod_cache[1] = log10(0.5 / err_0 * (err_2 / freq_b + err_1 / freq_a));
        }
        else
            lod_cache[0] = lod_cache[1] = 0;
        lod_cache[2] = log10((err_2 + minesp) / freq_b / (err_0 + minesp));
    }

    else if(archaic == '1'){
        if(selected){
            double err_3 = 3 - 2 * err_0;
            lod_cache[0] = -log10(freq_a * err_3);
            lod_cache[1] = -log10(2 * freq_a * freq_b * err_3);
            lod_cache[2] = -log10(freq_b * err_3);
        }
        else
            lod_cache[0] = lod_cache[1] = lod_cache[2] = 0;
    }

    else if(archaic == '2'){
        lod_cache[0] = log10((err_2 + minesp) / freq_a / (err_0 + minesp));
        if(selected){
            lod_cache[1] = log10(0.5 / err_0 * (err_1 / freq_b + err_2 / freq_a));
            lod_cache[2] = log10(err_1 / freq_b / err_0);
        }
        else
            lod_cache[1] = lod_cache[2] = 0;
    }

    else if(archaic == '9'){
        lod_cache[0] = 0;
        lod_cache[1] = 0;
        lod_cache[2] = 0;
    }
}

double Genotype_Reader::calculate_lod(char modern){
    if(modern == '9')
        return 0;
    return lod_cache[modern - '0'];
}

bool Genotype_Reader::yield_sample(char * &sample, int count){
    // iterate through samples
    // after execution, sample will point to the next sample
    // performed destructively
    // return true if sample points to a new sample

    // starting, initialize sample
    if(count == 0)
        sample = samples;
    // move to next sample
    else{
        for(; !is_break(sample); sample++);
        // recover tab in original
        if(count < num_samples)
            *sample = '\t';
        sample++;
    }
    if(count < num_samples){
        // modify sample to be null terminated (that sp points to)
        char * sp = sample;
        for(; !is_break(sp); sp++);
        *sp = '\0';
        return true;
    }
    return false;

}

bool is_break(const char * chr){
    return *chr == '\0' || *chr == '\t' || *chr == '\n';
}

int find_token(const char * query, const char * str){
    // find the token index of query in str
    // return -1 if not found
    int result = 0;
    if (str == nullptr)
        return -1;
    while(*str != '\0'){
        const char * qp, * sp;
        for(qp = query, sp = str; !is_break(sp) && !is_break(qp);
                qp++, sp++){
            if (*qp != *sp){
                break;
            }
        }
        //found
        if(is_break(qp) && is_break(sp)){
            return result;
        }

        // move str pointer to next token
        result++;
        for(; *sp != '\t' && *sp != '\0'; sp++);
        if(*sp == '\0')
            return -1;
        str = sp+1;
    }
    return -1;
}
