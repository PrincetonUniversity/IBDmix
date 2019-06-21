#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef unsigned long int ulnt;
int simpCharPar(const char c1, const char c2);

class VCF_File
{
    public:
        int chromosome;
        ulnt position;
        char reference;
        char alternative;
        int* genotypes;
        int number_individuals;
        FILE* input;
        char token[2048];

        VCF_File(FILE* in_file)
        {
            number_individuals = 0;
            input = in_file;
            // remove header
            while(fgetc(input) == '#'){
                if(fgetc(input) == '#')
                    purge_line();
                else{
                    char c;
                    while((c = fgetc(input)) != '\n')
                        if (c == '\t')
                            number_individuals++;
                    break;
                }
            }
            number_individuals -= 8; //remove other information
            genotypes = new int [number_individuals];
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

            fscanf(input, "%*s\t%*s\t%*s\t%*s\t");
            for(int i = 0; i < number_individuals; i++){
                fscanf(input, "%s", token);
                genotypes[i] = simpCharPar(token[0], token[2]);

                // convert anything above 2 to 9, or skip
                if(skip_non_informative && genotypes[i] == 9){
                    purge_line();
                    return false;
                }
            }
            return true;
        }

        ~VCF_File(){
            delete genotypes;
            fclose(input);
        }
};

int main(int argc, char *argv[])
{
    FILE * archaic_vcf, * modern_vcf, * output;
    int c;

    //read in arguments
    while ((c = getopt(argc, argv, "a:o:m:")) != -1) {
        switch (c) {
            case 'a': archaic_vcf = fopen(optarg, "r"); break;
            case 'm': modern_vcf = fopen(optarg, "r"); break;
            case 'o': output = fopen(optarg, "w"); break;
        }
    }

    VCF_File modern(modern_vcf);
    VCF_File archaic(archaic_vcf);

    bool recheck = false;
    //for each line in stdin
    while(modern.update()){

        //short getline with recheck for case when greater than
        while(recheck || archaic.update(true)){
            recheck = false;

            //less than position, copy all
            if(archaic.position < modern.position)
            {
                bool nonzero = false;
                for(int i = 0; i < archaic.number_individuals; i++)
                    if(archaic.genotypes[i] > 0){
                        nonzero = true;
                        break;
                    }

                if (nonzero){
                    fprintf(output, "%d\t%lu\t%c\t%c\t",
                            archaic.chromosome,
                            archaic.position,
                            archaic.reference,
                            archaic.alternative);
                    
                    for(int i = 0; i < archaic.number_individuals; i++)
                        fprintf(output, "%d\t", archaic.genotypes[i]);

                    //zero fill for stdin 
                    for(int i = 0; i < modern.number_individuals; i++)
                        fprintf(output, "0\t");

                    fprintf(output, "\n");
                }
            }
            //equal, have to check other conditions to write
            else if(archaic.position == modern.position)
            {
                if (archaic.reference == modern.reference &&
                        (archaic.alternative == '.' ||
                         archaic.alternative == modern.alternative))
                {
                    fprintf(output, "%d\t%lu\t%c\t%c\t",
                            archaic.chromosome,
                            archaic.position,
                            archaic.reference,
                            modern.alternative);
                    
                    for(int i = 0; i < archaic.number_individuals; i++)
                        fprintf(output, "%d\t", archaic.genotypes[i]);

                    for(int i = 0; i < modern.number_individuals; i++)
                        fprintf(output, "%d\t", modern.genotypes[i]);

                    fprintf(output, "\n");
                    break;
                }
                // advance modern to match legacy version
                break;
            }
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

inline int simpCharPar(const char c1, const char c2)
{
    // assumes genotype is either called or not, i.e. './0' is 9
    if(c1 == '0')
        return c2 == '0' ? 0 : 1;
    if(c1 == '1')
        return c2 == '0' ? 1 : 2;
    return 9;
}

