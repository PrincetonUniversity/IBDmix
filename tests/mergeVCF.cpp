#include <unistd.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

typedef unsigned long int ulnt;
// Compile: g++ -std=c++11 mergeVCF.cpp -o mergeVCF2

// Step 1. Filtering out multi-allelic SNVs and indels from archaic genomes
// Single Archaic genomes: $vcftools --gzvcf /net/akey/vol2/rcmccoy/archaic_vcf/Altai/chr1_mq25_mapab100.vcf.gz --out Altai_chr1 --remove-indels --min-alleles 1 --max-alleles 2 --recode
// Multiple Archaic genomes: $bcftools merge -m all --threads 2 ./Altai/chr1_mq25.vcf.gz ./Vindi/chr1_mq25.vcf.gz ./Denisova/chr1_mq25.vcf.gz | bcftools view --threads 2 --min-alleles 1 --max-alleles 2 -V indels -g ^miss -O z -o ./InterSectAVD_chr1_mq25.vcf.gz

// Step 2. Merge archaic vcf with vcf from 1000 genomes
// $gunzip InterSectAVD_chr1_mq25.vcf.gz
// $tabix -h ./ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 1:1-250000000 | ./mergeVCF -a ../Archaic/InterSectAVD_chr1_mq25.vcf -l 6260856 -d 3 -n 2504 -o InterSectAVD_1KGP3_chr1.txt

std::vector<std::string> split(const std::string& s, char seperator);
int simpCharPar(const char c);

int main(int argc, char *argv[]){
    ulnt position;
    int c;
    int refpos = 3;
    int altpos = 4;
    std::ifstream input;
    std::streambuf *buf;
    std::ofstream of;
    std::string inputLine, stdinLine;
    bool suppress = false;
    bool textInput = false;
    bool outputSet = false;
    
    //read in arguments
    while ((c = getopt(argc, argv, "a:o:s:")) != -1) {
        switch (c) {
            case 'a': input.open(optarg); break; // Archaic vcf from step 1
            case 'o':
                outputSet = true;
                of.open(optarg);
                buf = of.rdbuf();
                break;
            case 's': 
                suppress = (optarg[0] == 'y' || optarg[0] == 't');
                textInput = (optarg[0] == 't');
                break;
        }
    }

    //output to provided buffer or stdout
    if(outputSet == false)
        buf = std::cout.rdbuf();
    std::ostream output(buf);

    //remap the reference and alternate positions if text is provided
    if(textInput == true){
        refpos = 2;
        altpos = 3;
    }

    bool recheck = false;
    //for each line in stdin
    while(std::getline(std::cin, stdinLine)){
        //skip lines starting with #
        if(stdinLine[0] == '#')
            continue;

        std::vector<std::string> stdToks = split(stdinLine, '\t');
        
        //short getline with recheck for case when greater than
        while(recheck || std::getline(input, inputLine)){
            recheck = false;
            if(inputLine[0] == '#')
                continue;

            std::vector<std::string> inputToks = split(inputLine, '\t');

            // skip lines with multiple ref or alt nucleotides
            if(inputToks[refpos].length() > 1 || inputToks[altpos].length() > 1)
                continue;

            // skip lines with ./. in each entry
            if(!textInput){
                bool noneFound = true;
                for(int i = 9; i < inputToks.size(); i++){
                    int gt = simpCharPar(inputToks[i][0]) +
                        simpCharPar(inputToks[i][2]);
                    if(gt <= 2){ //valid entry
                        noneFound = false;
                        break;
                    }
                }

                if(noneFound)
                    continue;
            }

            //less than position, copy all, skip if suppressing
            if(stoi(inputToks[1]) < stoi(stdToks[1]))
            {
                //always suppress with text input!
                if(suppress == true)
                    continue;
                bool nonzero = false;
                for(int i = 9; i < inputToks.size(); i++)
                    if(simpCharPar(inputToks[i][0]) +
                            simpCharPar(inputToks[i][2]) > 0)
                    {
                        nonzero = true;
                        break;
                    }

                if(nonzero)
                {
                    output << inputToks[0] << "\t";
                    output << inputToks[1] << "\t";
                    output << inputToks[3] << "\t";
                    output << inputToks[4] << "\t";
                    
                    for(int i = 9; i < inputToks.size(); i++)
                        output << (simpCharPar(inputToks[i][0]) +
                                    simpCharPar(inputToks[i][2])) << "\t";

                    //zero fill for stdin 
                    for(int i = 9; i < stdToks.size(); i++)
                        output << "0\t";

                    output << "\n";
                                    
                }
            }
            //equal, have to check other conditions to write
            else if(stoi(inputToks[1]) == stoi(stdToks[1]))
            {
                if (stdToks[3].compare(inputToks[refpos]) == 0 &&
                        stdToks[4].length() == 1 &&
                        (inputToks[altpos][0] == '.' ||
                         inputToks[altpos].compare(stdToks[4]) == 0))
                {
                    output << inputToks[0] << "\t";
                    output << inputToks[1] << "\t";
                    output << inputToks[refpos] << "\t";
                    output << stdToks[4] << "\t";
                    
                    //note the -1 is due to the text output having an extra
                    //tab character at the end
                    if(textInput)
                        for(int i = 4; i < inputToks.size()-1; i++)
                            output << inputToks[i] << "\t";
                    else
                        for(int i = 9; i < inputToks.size(); i++)
                            output << (simpCharPar(inputToks[i][0]) +
                                        simpCharPar(inputToks[i][2])) << "\t";

                    //fill for stdin 
                    for(int i = 9; i < stdToks.size(); i++){
                        int gt = simpCharPar(stdToks[i][0]) +
                                    simpCharPar(stdToks[i][2]);
                        output << (gt <= 2 ? gt : 9) << "\t";
                    }

                    output << "\n";
                    break;
                }
                // advance VCF to match legacy version
                break;
            }
            else
            {
                recheck = true;
                break;
            }
        }
    }
    
    input.close();
    output.flush();
    if(outputSet == true)
        of.close();

    return 0;
}

std::vector<std::string> split(const std::string& s, char seperator)
{
    std::vector<std::string> output;
    std::string::size_type prev_pos = 0, pos = 0;

    while((pos = s.find(seperator, pos)) != std::string::npos)
    {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );
        output.push_back(substring);
        prev_pos = ++pos;
    }
    output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word

    return output;
}

int simpCharPar(const char c)
{
    switch(c)
    {
        case('1'): return 1;
        case('0'): return 0;
    }
    return 9;
}
