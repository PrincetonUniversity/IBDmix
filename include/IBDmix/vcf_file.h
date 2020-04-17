#pragma once

#include <iostream>
#include <sstream>
#include <string>

typedef unsigned long int ulnt;

class VCF_File
{
    private:
        std::istream* input;
        std::string buffer;
        std::istringstream iss;
    public:
        int chromosome;
        ulnt position;
        char reference;
        char alternative;
        std::string genotypes;
        std::string blank_line;
        int number_individuals;
        bool isvalid;

        VCF_File(std::istream* in_file, std::ostream &output);
        bool update(bool skip_non_informative=false);
        bool read_line(bool skip_non_informative=false);
};

