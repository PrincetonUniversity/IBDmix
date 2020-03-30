#include <string.h>
#include <getopt.h>
#include <iostream>
#include <fstream>

#include "IBDmix/vcf_file.h"

void print_options(void){
    std::cout << "--------------------------OPTIONS--------------------------\n"
    << "  -h, --help:       print help and exit\n"
    << "  -a, --archaic:    the archaic vcf file, uncompressed\n"
    << "  -m, --modern:     the modern vcf file, uncompressed\n"
    << "  -o, --output:     the merged, genotype file, uncompressed\n";
}

int main(int argc, char *argv[])
{
    FILE * archaic_vcf, * modern_vcf;
    std::string outfile = "-";
    int c;

    archaic_vcf = modern_vcf = nullptr;
    // read in arguments
    const option long_opts[] = {
        {"help", no_argument, nullptr, 'h'},
        {"archaic", required_argument, nullptr, 'a'},
        {"modern", required_argument, nullptr, 'm'},
        {"output", required_argument, nullptr, 'o'},
    };
    while ((c = getopt_long(argc, argv, "ha:o:m:", long_opts, nullptr)) != -1) {
        switch (c) {
            case 'h':
                std::cout << "Usage " << argv[0]
                    << " -a <archaic file> -m <modern file> -o <output file>\n";
                print_options();
                return 0;
            case 'a':
                archaic_vcf = strcmp(optarg, "-") == 0 ? stdin : fopen(optarg, "r");
                break;
            case 'm':
                modern_vcf = strcmp(optarg, "-") == 0 ? stdin : fopen(optarg, "r");
                break;
            case 'o':
                outfile = optarg;
                break;
        }
    }

    // error if any unset
    if (archaic_vcf == nullptr){
        fprintf(stderr,
                "Missing archaic vcf file input. Please provide valid '-a'\n");
        print_options();
        exit(1);
    }
    if (modern_vcf == nullptr){
        fprintf(stderr,
                "Missing modern vcf file input. Please provide valid '-m'\n");
        print_options();
        exit(1);
    }

    std::ofstream of;
    std::streambuf * buf;
    if (outfile == "-"){
        buf = std::cout.rdbuf();
    }
    else{
        of.open(outfile);
        buf = of.rdbuf();
    }
    std::ostream output(buf);

    // write header
    output << "chrom\tpos\tref\talt";
    // build files, this will print the headers as well (archaic first)
    VCF_File archaic(archaic_vcf, output);
    VCF_File modern(modern_vcf, output);
    // terminate header
    output << "\n";

    bool recheck = false;
    //for each line in modern
    while(modern.update()){

        //short update with recheck for case when archaic position is > modern
        while(recheck || archaic.update(true)){
            recheck = false;

            //less than position, copy archaic and use modern blank line
            if(archaic.position < modern.position)
            {
                // skip lines with no informative archaic GT
                bool nonzero = false;
                for(int i = 0; i < 2*archaic.number_individuals; i+=2)
                    if(archaic.genotypes[i] != '0'){
                        nonzero = true;
                        break;
                    }

                // found at least one informative archaic site
                if (nonzero){
                    // output archaic information and blank modern
                    output
                        << archaic.chromosome << '\t'
                        << archaic.position << '\t'
                        << archaic.reference << '\t'
                        << archaic.alternative << '\t';
                    output << archaic.genotypes;
                    output << modern.blank_line;
                    output << "\n";
                }
            }
            // equal, have to check other conditions to write
            else if(archaic.position == modern.position)
            {
                // check alleles are matching
                if (modern.isvalid && archaic.reference == modern.reference &&
                        (archaic.alternative == '.' ||
                         archaic.alternative == modern.alternative))
                {
                    output
                        << archaic.chromosome << '\t'
                        << archaic.position << '\t'
                        << archaic.reference << '\t'
                        << modern.alternative << '\t';
                    output << archaic.genotypes;
                    output << modern.genotypes;
                    output << "\n";
                }
                // advance both to match legacy version
                break;
            }
            // greater than, advance modern but keep archaic where it is
            else
            {
                recheck = true;
                break;
            }
        }
    }
    
    if(of.is_open())
        of.close();

    return 0;
}
