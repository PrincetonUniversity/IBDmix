#include <stdio.h>
#include <string.h>

#include <CLI/CLI.hpp>
#include <fstream>
#include <iostream>

#include "IBDmix/Genotype_Reader.h"
#include "IBDmix/IBD_Collection.h"
#include "IBDmix/IBD_Stack.h"

int main(int argc, char *argv[]) {
  CLI::App app{"Find probable IBD regions"};

  std::string genotype_file;
  app.add_option("-g,--genotype", genotype_file, "The genotype file")
      ->check(CLI::ExistingFile)
      ->required();

  std::string outfile = "-";
  app.add_option("-o,--output", outfile, "The output file location");

  std::string sample_file = "";
  app.add_option("-s,--sample", sample_file,
                 "File containing samples to "
                 "select from genotype.  Default to all samples in genotype.")
      ->check(CLI::ExistingFile);

  std::string archaic = "";
  app.add_option("-n,--archaic", archaic,
                 "Name of archaic sample, default"
                 " to first sample in genotype file");

  std::string mask_file = "";
  app.add_option("-r,--mask", mask_file,
                 "Mask of regions to 'remove'. "
                 "Regions in bed file have LOD set to 0")
      ->check(CLI::ExistingFile);

  int ma_threshold = 1;
  double LOD_threshold = 3.0, archaic_error = 0.01, modern_error_max = 0.0025,
         modern_error_prop = 2;
  // if output of should include additional stats about each region
  bool more_stats = false;
  bool inclusive_end = false;
  bool include_sites = false;
  bool include_lods = false;
  app.add_option("-d,--LOD-threshold", LOD_threshold,
                 "Threshold for emitting regions");
  app.add_option("-m,--minor-allele-count-threshold", ma_threshold,
                 "Threshold count for filtering minor alleles");
  app.add_option("-a,--archaic-error", archaic_error,
                 "Allele error rate for archaic DNA");
  app.add_option("-e,--modern-error-max", modern_error_max,
                 "Maximum allele error rate for modern samples");
  app.add_option("-c,--modern-error-proportion", modern_error_prop,
                 "Ratio between allele error rate and minor allele frequency");
  app.add_flag("-t,--more-stats", more_stats,
               "Flag to report additional region-level statistics");
  app.add_flag("-i,--inclusive-end", inclusive_end,
               "Change regions to be closed over [start, end]");
  app.add_flag("-w,--write-snps", include_sites,
               "Also include positions with positive LOD as a CSV list");
  app.add_flag("--write-lods", include_lods,
               "Also include LOD scores of positive LOD as a CSV list. "
               "Same order as SNPs.");

  CLI11_PARSE(app, argc, argv);

  std::ifstream genotype;
  genotype.open(genotype_file);

  std::ifstream sample;
  if (sample_file != "") sample.open(sample_file);
  std::ifstream mask;
  if (mask_file != "") mask.open(mask_file);

  // if output of end should be start, end) [exclusive end point]
  // or start, end] [inclusive end point; set exclusive to false]
  bool exclusive_end = !inclusive_end;

  std::ofstream of;
  std::streambuf *buf;
  if (outfile == "-") {
    buf = std::cout.rdbuf();
  } else {
    of.open(outfile);
    buf = of.rdbuf();
  }
  std::ostream output(buf);

  // write header
  output << "ID\tchrom\tstart\tend\tslod";

  Genotype_Reader reader(&genotype, &mask, archaic_error, modern_error_max,
                         modern_error_prop,
                         1e-200,  // minesp
                         ma_threshold);

  int num_samples = reader.initialize(sample, archaic);
  if (sample.is_open()) sample.close();

  IBD_Collection ibds(LOD_threshold, exclusive_end);

  ibds.initialize(reader);
  if (more_stats) ibds.add_recorder(IBD_Collection::Recorder::counts);
  if (include_sites) ibds.add_recorder(IBD_Collection::Recorder::sites);
  if (include_lods) ibds.add_recorder(IBD_Collection::Recorder::lods);

  ibds.writeHeader(output);
  output << '\n';

  while (reader.update()) ibds.update(reader, output);

  ibds.purge(output);

  genotype.close();
  if (mask.is_open()) mask.close();
  if (of.is_open()) of.close();
}
