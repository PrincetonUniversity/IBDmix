#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

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
        static const int buf_size = 1024;
        char buffer[buf_size];

        VCF_File(FILE* in_file, FILE* output)
        {
            // setup lines for subsequent reading, write individuals to output
            number_individuals = 0;
            input = in_file;
            char *individual_start = NULL;
            // remove header
            while(fgetc(input) == '#'){
                // header lines without needed information
                if(fgetc(input) == '#')
                    purge_line();
                // header with column names
                else{
                    bool remaining = true;
                    // if buffer hasn't read entire line
                    while(remaining){
                        // fill buffer
                        fgets(buffer, buf_size-1, input);
                        for(char *ptr = buffer; *ptr != '\0'; ptr++){
                            // count individuals, mark location for printing
                            if(*ptr == '\t'){
                                number_individuals++;
                                if(number_individuals > 8 &&
                                        individual_start == NULL)
                                    individual_start = ptr;
                            }
                            // line is read completely
                            if(*ptr == '\n'){
                                remaining = false;
                                *ptr = '\0';
                                break;
                            }
                        }
                        // write individuals to output
                        fprintf(output, "%s", individual_start);
                        individual_start = NULL;
                    }
                    break;
                }
            }
            number_individuals -= 8; //remove other columns
            // note: because we count tabs above, we subtract 8 (instead of 9)
            // as the last column has a newline instead of a tab

            // allocate output lines, fill in tabs
            genotypes = new char [number_individuals*2 + 1];
            blank_line = new char [number_individuals*2 + 1];
            for(int i = 1; i < 2*number_individuals; i+=2){
                genotypes[i] = '\t';
                blank_line[i] = '\t';
                blank_line[i-1] = '0';
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
            if(fscanf(input, "%i\t%lu\t%*s\t", &chromosome, &position) == EOF){
                chromosome = -1;
                return false;
            }

            // read in ref and alt alleles, skipping if > 1 character
            reference = fgetc(input);
            if (fgetc(input) != '\t'){
                purge_line();
                return false;
            }
            alternative = fgetc(input);
            if (fgetc(input) != '\t'){
                purge_line();
                return false;
            }

            // discard qual, filter, etc
            fscanf(input, "%*s\t%*s\t%*s\t%*s\t");
            bool remaining = true;
            // read in genotypes
            while(remaining){
                fgets(buffer, buf_size, input);
                char * ptr = buffer;
                int ind = 0;
                while(true){
                    genotypes[ind] = ptr[0] + ptr[2] - '0';
                    // ',' = '.' + '.' - '0'
                    genotypes[ind] = genotypes[ind] == ',' ? '9' : genotypes[ind];
                    if(skip_non_informative && genotypes[ind] == '9'){
                        if (strlen(buffer) == buf_size)
                            purge_line();
                        return false;
                    }
                    //read to next token
                    ind += 2;
                    while(*ptr != '\t' && *ptr != '\n' && *ptr != '\0'){
                        ptr++;
                    }
                    // end of line, break out of outer while loop
                    if(*ptr == '\n'){
                        remaining = false;
                        break;
                    }
                    // end of buffer, break out of inner while loop
                    if(*ptr == '\0')
                        break;
                    ptr++;
                }
            }
            return true;
        }

        ~VCF_File(){
            delete blank_line;
            delete genotypes;
            fclose(input);
        }
};

int main(int argc, char *argv[])
{
    FILE * archaic_vcf, * modern_vcf, * output;
    int c;

    archaic_vcf = modern_vcf = output = NULL;
    // read in arguments
    while ((c = getopt(argc, argv, "a:o:m:")) != -1) {
        switch (c) {
            case 'a': archaic_vcf = fopen(optarg, "r"); break;
            case 'm': modern_vcf = fopen(optarg, "r"); break;
            case 'o': output = fopen(optarg, "w"); break;
        }
    }

    // error if any unset
    if (archaic_vcf == NULL){
        printf("Missing archaic vcf file input. Please provide valid '-a'\n");
        exit(1);
    }
    if (output == NULL){
        printf("Missing output. Please provide '-o'\n");
        exit(1);
    }
    if (modern_vcf == NULL){
        printf("Missing modern vcf file input. Please provide valid '-m'\n");
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

            //less than position, copy archaic and use modern blankline
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
                if (archaic.reference == modern.reference &&
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
                    break;
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
