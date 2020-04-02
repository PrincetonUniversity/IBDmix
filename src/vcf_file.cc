#include "IBDmix/vcf_file.h"

VCF_File::VCF_File(FILE* in_file, std::ostream &output) {
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
            output << individual_start ;
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
    for(int i = 0; i < 2*number_individuals; i+=2){
        genotypes[i] = 'x';
        genotypes[i+1] = '\t';

        blank_line[i] = '0';
        blank_line[i+1] = '\t';
    }
    genotypes[number_individuals*2] = '\0';
    blank_line[number_individuals*2] = '\0';
}

void VCF_File::purge_line(void){
    fscanf(input, "%*[^\n]\n");
}

bool VCF_File::update(bool skip_non_informative){
    while(chromosome != -1 && !read_line(skip_non_informative))
        ;
    return chromosome != -1;
}

bool VCF_File::read_line(bool skip_non_informative){
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

VCF_File::~VCF_File(){
    delete []blank_line;
    delete []genotypes;
    if (buffer)
        free(buffer);
    fclose(input);
}
