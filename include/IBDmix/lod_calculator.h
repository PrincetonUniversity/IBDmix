#include <vector>
class LodCalculator {
 public:
  LodCalculator(double archaic_error = 0.01, double modern_error_max = 0.002,
                double modern_error_proportion = 2, double minesp = 1e-200)
      : archaic_error(archaic_error),
        modern_error_max(modern_error_max),
        modern_error_proportion(modern_error_proportion),
        minesp(minesp),
        lod_cache(3) {}

  double get_modern_error(double frequency) const;
  void update_lod_cache(char archaic, double freq_b, bool selected = true);
  double calculate_lod(char modern) const;

 private:
  double archaic_error;
  double modern_error_max;
  double modern_error_proportion;
  double minesp;
  std::vector<double> lod_cache;
};
