#include <getopt.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "IBD_Collection.hpp"
#include "Genotype_Reader.hpp"
#include "IBD_Stack.hpp"

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
                    << " -c <modern error proportion>\n";
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
            reader);

    while(reader.update())
        ibds.update(reader, output);
    //TODO write tests, acceptance on real data

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
