#pragma once
#include <stdio.h>

int find_token(const char * query, const char * str);

class Genotype_Reader{
    private:
        FILE * genotype;
        FILE * mask;
        int *sample_to_index, num_samples, archaic_index, minor_allele_cutoff,
            mask_chromosome;
        unsigned long int mask_start, mask_end;
        char *buffer;
        char *samples;
        size_t buf_size;
        double archaic_error, modern_error_max, modern_error_proportion,
               minesp, *lod_cache;

        void determine_sample_mapping(const char * sample_line);
        void find_archaic(const char * archaic,
                const char * sample_line);
        bool in_mask(void);
        void process_line_buffer(bool selected);
        bool get_frequency(double &frequency);
        double get_modern_error(double frequency);
        double calculate_lod(char modern);
        void update_lod_cache(char archaic, double freq_b, double modern_error,
                bool selected=true);
    public:
        Genotype_Reader(FILE * genotype, FILE * mask=nullptr,
                double archaic_error=0.01, double modern_error_max=0.002,
                double modern_error_proportion=2, double minesp=1e-200,
                int minor_allele_cutoff=1);
        ~Genotype_Reader();

        int chromosome;
        unsigned long int position;
        double *lod_scores;

        int initialize(FILE * samples=nullptr,
                const char * archaic=nullptr);
        bool update(void);
        bool yield_sample(char * &sample, int count);

};
