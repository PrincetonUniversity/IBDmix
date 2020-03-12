#include <getopt.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "IBDmix/IBD_Collection.h"
#include "IBDmix/Genotype_Reader.h"
#include "IBDmix/IBD_Stack.h"

void print_options(void){
    std::cout << "--------------------------OPTIONS--------------------------\n"
        << "  -h, --help:       print help and exit\n"                // <- max
        << "\n"
        << "  -g, --genotype:   input genotype file for one chromosome\n"
        << "                      REQUIRED\n"
        << "\n"
        << "  -o, --output:     output file\n"
        << "                      REQUIRED\n"
        << "\n"
        << "  -s, --sample:     sample file, one sample name per line\n"
        << "                      DEFAULT: use all samples\n"
        << "\n"
        << "  -n, --archaic:    name of archaic sample\n"
        << "                      DEFAULT: first sample taken as archaic\n"
        << "\n"
        << "  -r, --mask:       mask file, sorted, merged bed file\n"
        << "                      DEFAULT: use all positions\n"
        << "\n"
        << "  -d, --LOD-threshold:\n"
        << "                    threshold of log(odds) to emit into \n"
        << "                    output file\n"
        << "                      DEFAULT: 3.0\n"
        << "\n"
        << "  -m, --minor-allele-count-threshold:\n"
        << "                    threshold for filtering minor alleles\n"
        << "                      DEFAULT: 1\n"
        << "\n"
        << "  -a, --archaic-error:\n"
        << "                    allele error rate for archaic DNA\n"
        << "                      DEFAULT: 0.01\n"
        << "\n"
        << "  -e, --modern-error-max\n"
        << "                    maximum allele error rate for modern humans\n"
        << "                      DEFAULT: 0.0025\n"
        << "\n"
        << "  -c, --modern-error-proportion\n"
        << "                    ratio between allele error rate and minor\n"
        << "                    allele frequency for modern humans\n"
        << "                      DEFAULT: 2\n"
        << "\n"
        << "  -i, --inclusive-end\n"
        << "                    switch to report the end position as the\n"
        << "                    first position of decreasing LOD score.\n"
        << "                      DEFAULT: end is *next* position in genotype file\n"
        << "\n"
        << "  -t, --more-stats\n"
        << "                    report additional columns detailing execution\n"
        << "                    sites: total number of sites in genotype file in region\n"
        << "                    mask and maf: sites which are both in the masked region and fail MAF cutoff\n"
        << "                    in_mask: sites in the mask but pass MAF\n"
        << "                    maf_low: sites failing MAF with low frequency\n"
        << "                    maf_high: sites failing MAF with high frequency\n"
        << "                    rec_2_0: sites failing mask or MAF but still considered due to a genotype of 2 in archaic and 0 in modern\n"
        << "                    rec_0_2: same as rec_2_0 but with 0 in archaic and 2 in modern\n"
        << "                      DEFAULT: off\n"
        ;
}

int main(int argc, char *argv[]){
    FILE *genotype = nullptr,
         *sample = nullptr,
         *mask = nullptr,
         *output = nullptr;
    int ma_threshold = 1;
    char * archaic = nullptr;
    double LOD_threshold = 3.0,
           archaic_error = 0.01,
           modern_error_max = 0.0025,
           modern_error_prop = 2;
    // if output of end should be start, end) [exclusive end point]
    // or start, end] [inclusive end point; set exclusive to false]
    bool exclusive_end = true;
    // if output of should include additional stats about each region
    bool more_stats = false;

    static const struct option long_opts[] = {
        {"help", no_argument, nullptr, 'h'},
        {"genotype", required_argument, nullptr, 'g'},
        {"output", required_argument, nullptr, 'o'},
        {"sample", required_argument, nullptr, 's'},
        {"archaic", required_argument, nullptr, 'n'},
        {"mask-file", required_argument, nullptr, 'r'},
        {"LOD-threshold", required_argument, nullptr, 'd'},
        {"minor-allele-count-threshold", required_argument, nullptr, 'm'},
        {"archaic-error", required_argument, nullptr, 'a'},
        {"modern-error-max", required_argument, nullptr, 'e'},
        {"modern-error-proportion", required_argument, nullptr, 'c'},
        {"inclusive-end", no_argument, nullptr, 'i'},
        {"more-stats", no_argument, nullptr, 't'},
        {0, 0, 0, 0}
    };
    char c;
    while ((c = getopt_long(argc, argv, "hg:s:a:o:d:m:a:e:c:r:",
                    long_opts, nullptr)) != -1) {
        switch (c) {
            case 'h':
                std::cout << "Usage " << argv[0]
                    << " -g <genotype file> -s <sample file> -r <mask file>\n"
                    << " -o <output> -d <LOD threshold> -m <MA cutoff>\n"
                    << " -a <archaic error> -e <modern error max>\n"
                    << " -c <modern error proportion> [-i] [-t]\n";
                print_options();
                return 0;
            case 'g': 
                genotype = strcmp(optarg, "-") == 0 ? stdin : fopen(optarg, "r");
                break;
            case 'o':
                output = strcmp(optarg, "-") == 0 ? stdin : fopen(optarg, "w");
                break;
            case 's':
                sample = strcmp(optarg, "-") == 0 ? stdout : fopen(optarg, "r");
                break;
            case 'n':
                archaic = new char[strlen(optarg)];
                strcpy(archaic, optarg);
                break;
            case 'r': mask = fopen(optarg, "r"); break;
            case 'd': LOD_threshold = atof(optarg); break;
            case 'm': ma_threshold = atoi(optarg); break;
            case 'a': archaic_error = atof(optarg); break;
            case 'e': modern_error_max = atof(optarg); break;
            case 'c': modern_error_prop = atof(optarg); break;
            case 'i': exclusive_end = false; break;
            case 't': more_stats = true; break;
        }
    }

    if (genotype == nullptr){
        printf("Missing genotype file input. Please provide valid '-g'\n");
        print_options();
        exit(1);
    }
    if (output == nullptr){
        printf("Missing output file. Please provide valid '-o'\n");
        print_options();
        exit(1);
    }

    // write header
    fprintf(output, "%s\t%s\t%s\t%s\t%s",
            "ID", "chrom", "start", "end", "slod");
    if (more_stats == true)
        fprintf(output, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                "sites", "mask and maf", "in_mask", "maf_low",
                "maf_high", "rec_2_0", "rec_0_2");
    fprintf(output, "\n");
    Genotype_Reader reader = Genotype_Reader(
            genotype,
            mask,
            archaic_error,
            modern_error_max,
            modern_error_prop,
            1e-200,  // minesp
            ma_threshold);
    int num_samples = reader.initialize(sample, archaic);
    if(sample != nullptr)
        fclose(sample);

    IBD_Collection ibds;
    ibds.initialize(
            num_samples,
            LOD_threshold,
            reader,
            exclusive_end,
            more_stats);

    while(reader.update())
        ibds.update(reader, output);

    ibds.purge(output);

    if (archaic != nullptr)
        free(archaic);
    if(genotype != nullptr)
        fclose(genotype);
    if(mask != nullptr)
        fclose(mask);
    if(output != nullptr)
        fclose(output);

    free_stack();
}
