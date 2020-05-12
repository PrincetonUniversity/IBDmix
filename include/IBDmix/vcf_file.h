#pragma once

#include <iostream>
#include <sstream>
#include <string>

class VCF_File {
 private:
  std::istream *input;
  std::string buffer;
  std::istringstream iss;

  bool simpleParse(const char *start);
  bool complexParse(const char *start, int gtInd);
  bool parse(const char *start, std::string format);

 public:
  int chromosome;
  uint64_t position;
  char reference;
  char alternative;
  std::string genotypes;
  std::string blank_line;
  int number_individuals;
  bool isvalid;

  VCF_File(std::istream *in_file, std::ostream &output);
  bool update(bool skip_non_informative = false);
  bool read_line(bool skip_non_informative = false);
};
