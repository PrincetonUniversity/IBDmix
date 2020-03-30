#pragma once

#include <stdio.h>
#include <stdlib.h>
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

        VCF_File(FILE* in_file, std::ostream &output);
        void purge_line();
        bool update(bool skip_non_informative=false);
        bool read_line(bool skip_non_informative=false);
        ~VCF_File();
};

