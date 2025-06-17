#pragma once

#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>

class VCF_File {
 public:
  VCF_File(std::istream *in_file, std::ostream &output);
  bool update(bool skip_non_informative = false);
  bool read_line(bool skip_non_informative = false);

  const std::string &getBlank() const { return blank_line; }
  const std::string &getGenotypes() const { return genotypes; }

  const std::string &getChromosome() const { return chromosome; }
  uint64_t getPosition() const { return position; }
  char getReference() const { return reference; }
  char getAlternative() const { return alternative; }
  int getCount() const { return number_individuals; }
  bool isValid() const { return isvalid; }

 private:
  std::istream *input;
  std::string buffer;
  std::istringstream iss;

  std::string genotypes;
  std::string blank_line;

  std::string chromosome;
  uint64_t position;
  char reference;
  char alternative;
  int number_individuals;
  bool isvalid;

  bool simpleParse(const char *start);
  bool complexParse(const char *start, int gtInd);
  bool parse(const char *start, std::string format);
};
