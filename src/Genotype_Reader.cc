#include "IBDmix/Genotype_Reader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>

int Genotype_Reader::initialize(std::istream &samples, std::string archaic) {
  // using samples list and header line, determine number of samples
  // and mapping from position to sample number
  std::getline(*genotype, buffer);
  std::istringstream iss(buffer);
  int result = sample_mapper.initialize(iss, samples, archaic);

  lod_scores.resize(result);
  recover_type.resize(result);
  return result;
}

bool Genotype_Reader::update() {
  // read next line of input file
  // update the lod_scores for reading, handling masks
  std::getline(*genotype, buffer);
  iss.clear();
  iss.str(buffer);

  // return false if the file is read fully
  if (!(iss >> chromosome && iss >> position)) return false;

  iss >> token;  // ref
  ref = token[0];
  iss >> token;  // alt
  alt = token[0];
  line_filtering = 0;
  // selected indicates if the line should have its lod calculated
  // set to false if one of the following occurs:
  // - in a masked region
  // - fails to meet allele cutoff
  // If selected is false, lod = 0, unless archaic = (0,2) and modern = (2,0)
  bool selected = !mask.in_mask(chromosome, position);
  if (!selected) line_filtering |= IN_MASK;

  // find the 4th tab and erase from buffer
  std::string::size_type ind = buffer.find('\t');
  for (int i = 1; i < 4; ++i) ind = buffer.find('\t', ind + 1);
  buffer.erase(0, ind + 1);
  process_line_buffer(selected);
  return true;
}

void Genotype_Reader::process_line_buffer(bool selected) {
  // assume buffer is loaded with tab-separated character in {0, 1, 2, 9}
  // from a genotype file.  Using sample_mapper, fill in
  // the lod_scores array with appropriate values

  // throughout, *2 to skip tabs
  archaic = buffer[sample_mapper.archaic_index * 2];
  selected &= find_frequency();

  calculator.update_lod_cache(archaic, allele_frequency, selected);

  for (int i = 0; i < sample_mapper.size(); i++) {
    lod_scores[i] =
        calculator.calculate_lod(buffer[sample_mapper.sample_to_index[i] * 2]);
    recover_type[i] = 0;
  }

  // udpate recover type
  if (!selected && archaic == '0') {
    for (int i = 0; i < sample_mapper.size(); i++)
      if (buffer[sample_mapper.sample_to_index[i] * 2] == '2')
        recover_type[i] = RECOVER_0_2;
  } else if (!selected && archaic == '2') {
    for (int i = 0; i < sample_mapper.size(); i++)
      if (buffer[sample_mapper.sample_to_index[i] * 2] == '0')
        recover_type[i] = RECOVER_2_0;
  }
}

bool Genotype_Reader::find_frequency() {
  // determine the observed frequency of alternative alleles
  // Returns true if enough counts were observed above the cutoff value
  int total_counts = 0;
  double alt_counts = 0;
  char current;
  bool select = true;
  for (int i = 0; i < sample_mapper.size(); i++)
    if ((current = buffer[sample_mapper.sample_to_index[i] * 2]) != '9') {
      total_counts += 2;
      alt_counts += current - '0';
    }

  if (alt_counts <= minor_allele_cutoff) {  // not enough counts
    select = false;
    line_filtering |= MAF_LOW;
  }
  if (total_counts - alt_counts <= minor_allele_cutoff) {  // too many
    select = false;
    line_filtering |= MAF_HIGH;
  }

  if (total_counts == 0)
    allele_frequency = 0;
  else
    allele_frequency = alt_counts / total_counts;
  return select;
}

const std::vector<std::string> &Genotype_Reader::get_samples() const {
  return sample_mapper.samples;
}
