#pragma once
#include <iostream>
#include <sstream>
#include <string>

class Mask_Reader{
    std::string chromosome = "";
    unsigned long start, end;
    std::istream *mask = nullptr;
    void readline();

public:
    Mask_Reader(std::istream *mask) : mask(mask) {readline();};
    bool in_mask(const std::string &chrom, unsigned long position);
};
