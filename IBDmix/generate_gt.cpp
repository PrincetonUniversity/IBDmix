#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>

typedef unsigned long int ulnt;

class VCF_File
{
    public:
        int chromosome;
        ulnt position;
        char reference;
        char alternative;
        char* genotypes;
        char* blank_line;
        int number_individuals;
        FILE* input;
        char* buffer;
        size_t len;
        bool isvalid;

        VCF_File(FILE* in_file, FILE* output)
        {
            // setup lines for subsequent reading, write individuals to output
            number_individuals = 0;
            input = in_file;
            char *individual_start = NULL;
            buffer = NULL;
            len = 0;
            char *ptr;
            chromosome = 0;
            // remove header
            while(fgetc(input) == '#'){
                // header lines without needed information
                if(fgetc(input) == '#')
                    purge_line();
                // header with column names
                else{
                    getline(&buffer, &len, input);
                    for(ptr = buffer; *ptr != '\0'; ptr++){
                        // count individuals, mark location for printing
                        if(*ptr == '\t'){
                            number_individuals++;
                            if(individual_start == NULL &&
                                    number_individuals > 8)
                                individual_start = ptr;
                        }
                        if(*ptr == '\n'){
                            *ptr = '\0';
                            // handle tab at line end for header and remove an indiv
                            if(*(ptr-1) == '\t'){
                                *(ptr-1) = '\0';
                                --number_individuals;
                            }
                            break;
                        }
                    }
                    // write individuals to output
                    fprintf(output, "%s", individual_start);
                    break;
                }
            }
            if(number_individuals <= 8){
                fprintf(stderr, "Ill-formed file, unable to parse header\n");
                exit(1);
            }

            number_individuals -= 8; //remove other columns
            // note: because we count tabs above, we subtract 8 (instead of 9)
            // as the last column has a newline instead of a tab

            // allocate output lines, fill in tabs
            genotypes = new char [number_individuals*2 + 1];
            blank_line = new char [number_individuals*2 + 1];
            for(int i = 0; i <= 2*number_individuals; i+=2){
                genotypes[i] = 'x';
                genotypes[i+1] = '\t';

                blank_line[i] = '0';
                blank_line[i+1] = '\t';
            }
            genotypes[number_individuals*2] = '\0';
            blank_line[number_individuals*2] = '\0';
        }

        void purge_line(void){
            fscanf(input, "%*[^\n]\n");
        }

        bool update(bool skip_non_informative=false){
            while(chromosome != -1 && !read_line(skip_non_informative))
                ;
            return chromosome != -1;
        }

        bool read_line(bool skip_non_informative=false){
            isvalid = true;
            if(fscanf(input, "%i\t%lu\t%*s\t", &chromosome, &position) == EOF){
                chromosome = -1;
                return false;
            }

            // read in ref and alt alleles, skipping if > 1 character
            reference = fgetc(input);
            if (fgetc(input) != '\t'){
                purge_line();
                isvalid = false;
                return !skip_non_informative;  // this will allow checks if not skipping
            }
            alternative = fgetc(input);
            if (fgetc(input) != '\t'){
                purge_line();
                isvalid = false;
                return !skip_non_informative;
            }

            // discard qual, filter, etc
            fscanf(input, "%*s\t%*s\t%*s\t%*s\t");
            // read in genotypes
            getline(&buffer, &len, input);
            int ind = 0;
            char * ptr = buffer;
            bool none_valid = true;
            while(true){
                genotypes[ind] = ptr[0] + ptr[2] - '0';
                // ',' = '.' + '.' - '0'
                genotypes[ind] = genotypes[ind] == ',' ? '9' : genotypes[ind];
                if(none_valid && genotypes[ind] != '9')
                    none_valid = false;
                // read to next token
                ind += 2;
                while(*ptr != '\t' && *ptr != '\0'){
                    ptr++;
                }
                // the check for newlines will handle tabs before newlines
                if(*ptr == '\0' || *(++ptr) == '\n')
                    break;
            }
            // return true of skip is false or
            // if skip is false but at least one informative found
            return !skip_non_informative || !none_valid;
        }

        ~VCF_File(){
            delete []blank_line;
            delete []genotypes;
            if (buffer)
                free(buffer);
            fclose(input);
        }
};

void print_options(void){
    std::cout << "--------------------------OPTIONS--------------------------\n"
    << "  -h, --help:       print help and exit\n"
    << "  -a, --archaic:    the archaic vcf file, uncompressed\n"
    << "  -m, --modern:     the modern vcf file, uncompressed\n"
    << "  -o, --output:     the merged, genotype file, uncompressed\n";
}

int main(int argc, char *argv[])
{
    FILE * archaic_vcf, * modern_vcf, * output;
    int c;

    archaic_vcf = modern_vcf = output = NULL;
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
                output = strcmp(optarg, "-") == 0 ? stdout : fopen(optarg, "w");
                break;
        }
    }

    // error if any unset
    if (archaic_vcf == NULL){
        fprintf(stderr,
                "Missing archaic vcf file input. Please provide valid '-a'\n");
        print_options();
        exit(1);
    }
    if (output == NULL){
        fprintf(stderr,
                "Missing output. Please provide '-o'\n");
        print_options();
        exit(1);
    }
    if (modern_vcf == NULL){
        fprintf(stderr,
                "Missing modern vcf file input. Please provide valid '-m'\n");
        print_options();
        exit(1);
    }

    // write header
    fprintf(output, "chrom\tpos\tref\talt");
    // build files, this will print the headers as well (archaic first)
    VCF_File archaic(archaic_vcf, output);
    VCF_File modern(modern_vcf, output);
    // terminate header
    fprintf(output, "\n");

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
                    fprintf(output, "%d\t%lu\t%c\t%c\t",
                            archaic.chromosome,
                            archaic.position,
                            archaic.reference,
                            archaic.alternative);
                    fprintf(output, "%s", archaic.genotypes);
                    fprintf(output, "%s", modern.blank_line);
                    fprintf(output, "\n");
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
                    fprintf(output, "%d\t%lu\t%c\t%c\t",
                            archaic.chromosome,
                            archaic.position,
                            archaic.reference,
                            modern.alternative);
                    fprintf(output, "%s", archaic.genotypes);
                    fprintf(output, "%s", modern.genotypes);
                    fprintf(output, "\n");
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
    
    fclose(output);

    return 0;
}
