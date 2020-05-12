#include "IBDmix/Genotype_Reader.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>

Genotype_Reader::Genotype_Reader(std::istream *genotype, std::istream *mask,
                                 double archaic_error, double modern_error_max,
                                 double modern_error_proportion, double minesp,
                                 int minor_allele_cutoff)
    : genotype(genotype),
      mask(mask),
      archaic_error(archaic_error),
      modern_error_max(modern_error_max),
      modern_error_proportion(modern_error_proportion),
      minesp(minesp),
      minor_allele_cutoff(minor_allele_cutoff) {
  lod_cache.resize(3);
}

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
  bool selected = !mask.in_mask(std::to_string(chromosome), position);
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
  double modern_error = get_modern_error();

  update_lod_cache(archaic, allele_frequency, modern_error, selected);

  for (int i = 0; i < sample_mapper.size(); i++) {
    lod_scores[i] = calculate_lod(buffer[sample_mapper.sample_to_index[i] * 2]);
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
  int total_counts = 0, alt_counts = 0;
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
    allele_frequency = static_cast<double>(alt_counts) / total_counts;
  return select;
}

double Genotype_Reader::get_modern_error() {
  double frequency = allele_frequency;
  if (frequency > 0.5)  // convert to minor frequency
    frequency = 1 - frequency;
  double prop_error = frequency * modern_error_proportion;

  return prop_error < modern_error_max ? prop_error : modern_error_max;
}

void Genotype_Reader::update_lod_cache(char archaic, double freq_b,
                                       double modern_error, bool selected) {
  // update the lod cache array based on values along genotype file line
  // TODO(troycomi) should minesp be added to more terms? (e.g. 2,2)
  double freq_a = 1 - freq_b;
  double err_0 = 1 - archaic_error * (1 - archaic_error);
  double err_1 =
      (1 - archaic_error) * (1 - modern_error) + archaic_error * modern_error;
  double err_2 =
      (1 - archaic_error) * modern_error + archaic_error * (1 - modern_error);

  if (archaic == '0') {
    if (selected) {
      lod_cache[0] = log10(err_1 / freq_a / err_0);
      lod_cache[1] = log10(0.5 / err_0 * (err_2 / freq_b + err_1 / freq_a));
    } else {
      lod_cache[0] = lod_cache[1] = 0;
    }
    lod_cache[2] = log10((err_2 + minesp) / freq_b / (err_0 + minesp));
  } else if (archaic == '1') {
    if (selected) {
      double err_3 = 3 - 2 * err_0;
      lod_cache[0] = -log10(freq_a * err_3);
      lod_cache[1] = -log10(2 * freq_a * freq_b * err_3);
      lod_cache[2] = -log10(freq_b * err_3);
    } else {
      lod_cache[0] = lod_cache[1] = lod_cache[2] = 0;
    }
  } else if (archaic == '2') {
    lod_cache[0] = log10((err_2 + minesp) / freq_a / (err_0 + minesp));
    if (selected) {
      lod_cache[1] = log10(0.5 / err_0 * (err_1 / freq_b + err_2 / freq_a));
      lod_cache[2] = log10(err_1 / freq_b / err_0);
    } else {
      lod_cache[1] = lod_cache[2] = 0;
    }
  } else if (archaic == '9') {
    lod_cache[0] = lod_cache[1] = lod_cache[2] = 0;
  }
}

double Genotype_Reader::calculate_lod(char modern) {
  if (modern == '9') return 0;
  return lod_cache[modern - '0'];
}

const std::vector<std::string> &Genotype_Reader::get_samples() const {
  return sample_mapper.samples;
}
