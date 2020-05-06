#ifndef INCLUDE_IBDMIX_MASK_READER_H_
#define INCLUDE_IBDMIX_MASK_READER_H_

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

class Mask_Reader{
 private:
    std::string chromosome = "";
    uint64_t start, end;
    std::istream *mask = nullptr;
    void readline();

 public:
    explicit Mask_Reader(std::istream *mask) : mask(mask) {readline();}
    bool in_mask(const std::string &chrom, uint64_t position);
};

#endif  // INCLUDE_IBDMIX_MASK_READER_H_
