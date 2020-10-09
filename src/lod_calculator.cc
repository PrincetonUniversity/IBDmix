#include "IBDmix/lod_calculator.h"

#include <math.h>

double LodCalculator::get_modern_error(double frequency) const {
  if (frequency > 0.5)  // convert to minor frequency
    frequency = 1 - frequency;
  double prop_error = frequency * modern_error_proportion;

  return prop_error < modern_error_max ? prop_error : modern_error_max;
}

void LodCalculator::update_lod_cache(char archaic, double freq_b,
                                     bool selected) {
  // update the lod cache array based on values along genotype file line
  // TODO(troycomi) should minesp be added to more terms? (e.g. 2,2)
  double modern_error = get_modern_error(freq_b);
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

double LodCalculator::calculate_lod(char modern) const {
  if (modern == '9') return 0;
  return lod_cache[modern - '0'];
}
