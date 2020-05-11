#include "IBDmix/vcf_file.h"

#include <string.h>

#include <stdexcept>

VCF_File::VCF_File(std::istream* in_file, std::ostream &output) :
    input(in_file) {
    // setup lines for subsequent reading, write individuals to output
    number_individuals = 0;
    chromosome = 0;
    // remove header
    while (std::getline(*input, buffer)) {
        // header lines without needed information
        if (buffer[1] == '#') {
            continue;
        // header with column names
        } else if (buffer[0] == '#') {
            iss.str(buffer);
            std::string token;
            // read in first 9 columns (not needed)
            iss >> token;  // CHROM
            iss >> token;  // POS
            iss >> token;  // ID
            iss >> token;  // REF
            iss >> token;  // ALT
            iss >> token;  // QUAL
            iss >> token;  // FILTER
            iss >> token;  // INFO
            iss >> token;  // FORMAT
            // write indivs prepend with \t, count number
            while (iss >> token) {
                ++number_individuals;
                output << '\t' << token;
            }
            break;
        }
    }

    if (number_individuals == 0) {
        std::cerr << "Ill-formed file, unable to parse header\n";
        exit(1);
    }

    // allocate output lines, fill in tabs
    genotypes.resize(number_individuals*2);
    blank_line.resize(number_individuals*2);
    for (int i = 0; i < 2*number_individuals; i+=2) {
        genotypes[i] = 'x';
        genotypes[i+1] = '\t';

        blank_line[i] = '0';
        blank_line[i+1] = '\t';
    }
}

bool VCF_File::update(bool skip_non_informative) {
    while (chromosome != -1 && !read_line(skip_non_informative)) {}
    return chromosome != -1;
}

bool VCF_File::read_line(bool skip_non_informative) {
    isvalid = true;
    if (!std::getline(*input, buffer)) {
        chromosome = -1;
        return false;
    }
    iss.clear();
    iss.str(buffer);
    if (!(iss >> chromosome && iss >> position))
        return false;

    // read in ref and alt alleles, skipping if > 1 character
    std::string token;
    iss >> token;  // ID
    iss >> token;  // REF
    if (token.size() != 1) {
        isvalid = false;
        return !skip_non_informative;  // this will allow checks if not skipping
    }
    reference = token[0];
    iss >> token;  // ALT
    if (token.size() != 1) {
        isvalid = false;
        return !skip_non_informative;  // this will allow checks if not skipping
    }
    alternative = token[0];

    // discard qual, filter, etc
    iss >> token;  // QUAL
    iss >> token;  // FILTER
    iss >> token;  // INFO
    iss >> token;  // FORMAT

    // skip to sample
    std::string::size_type ind = buffer.find('\t');
    for (int i = 1; i < 9; ++i)
        ind = buffer.find('\t', ind+1);

    bool none_valid = parse(buffer.c_str() + ind + 1, token);

    // return true if skip is false or
    // if skip is false but at least one informative found
    return !skip_non_informative || !none_valid;
}

bool VCF_File::parse(const char *start, std::string format) {
    // check if format is GT, otherwise need to parse more carefully
    if (format == "GT") {
        return simpleParse(start);
    } else {
        // format is complex, split by : finding index of GT
        int ind = 0;
        const char *fmt = format.c_str();
        for (;;) {
            if (strncmp(fmt, "GT", 2) == 0)
                return complexParse(start, ind);
            ++ind;  // not found, onto next fmt
            while (*fmt != ':' && *fmt != '\0')
                ++fmt;
            if (*fmt == '\0') {
                throw std::invalid_argument("FORMAT must contain GT");
            }
            ++fmt;
        }
    }
}

bool VCF_File::simpleParse(const char *start) {
    bool none_valid = true;
    unsigned int ind = 0;
    while (ind < genotypes.size()) {
        genotypes[ind] = start[0] + start[2] - '0';
        // ',' = '.' + '.' - '0'
        genotypes[ind] = genotypes[ind] == ',' ? '9' : genotypes[ind];
        if (none_valid && genotypes[ind] != '9')
            none_valid = false;
        ind += 2;
        start += 4;  // gt (2) separator and tab
    }
    return none_valid;
}

bool VCF_File::complexParse(const char *start, int gtInd) {
    bool none_valid = true;
    unsigned int ind = 0;
    while (ind < genotypes.size()) {
        // move to the gtInd'th :
        for (int i = 0; i < gtInd; ++i) {
            for (; *start != ':'; ++start) {}
            ++start;
        }

        genotypes[ind] = start[0] + start[2] - '0';
        // ',' = '.' + '.' - '0'
        genotypes[ind] = genotypes[ind] == ',' ? '9' : genotypes[ind];
        if (none_valid && genotypes[ind] != '9')
            none_valid = false;

        ind += 2;

        // move to next tab or null
        while (*start != '\t' && *start != '\0')
            ++start;
        ++start;
    }
    return none_valid;
}
