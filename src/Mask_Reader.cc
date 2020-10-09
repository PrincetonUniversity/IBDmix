#include "IBDmix/Mask_Reader.h"

bool Mask_Reader::in_mask(const std::string &chrom, uint64_t position) {
  // test if chrom/position is in mask file
  // assume queries are sorted in same order as mask file!
  // true if position is in (start, end]
  if (mask == nullptr) return false;

  for (;;) {
    if (chrom == chromosome) {
      if (position <= start)
        return false;
      else if (end < position)
        readline();
      else
        return true;  // start < position <= end
    } else if (chromosome == "") {
      return false;
      // this assumes same order. may fail with numeric vs lexigraphic sort
    } else if (chromosome < chrom) {
      readline();
    } else {
      return false;
    }
  }
}

void Mask_Reader::readline() {
  if (mask == nullptr) return;
  std::string line;
  if (std::getline(*mask, line)) {
    std::istringstream iss(line);
    if (!(iss >> chromosome >> start >> end)) {
      chromosome = "";
      throw std::invalid_argument("Unable to read mask file " + line);
    }
  } else {
    chromosome = "";
  }
}
