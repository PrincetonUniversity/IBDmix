#include <CLI/CLI.hpp>
#include <fstream>
#include <iostream>

#include "IBDmix/vcf_file.h"

int main(int argc, char *argv[]) {
  CLI::App app{"Produce genotype files from vcfs"};

  std::string archaic_file;
  app.add_option("-a,--archaic", archaic_file, "The archaic sample vcf")
      ->check(CLI::ExistingFile)
      ->required();

  std::string modern_file;
  app.add_option("-m,--modern", modern_file, "The modern sample vcf")
      ->check(CLI::ExistingFile)
      ->required();

  std::string outfile = "-";
  app.add_option("-o,--output", outfile, "The output file location");

  CLI11_PARSE(app, argc, argv);

  std::ofstream of;
  std::streambuf *buf;
  if (outfile == "-") {
    buf = std::cout.rdbuf();
  } else {
    of.open(outfile);
    buf = of.rdbuf();
  }
  std::ostream output(buf);

  std::ifstream archaic_vcf, modern_vcf;
  archaic_vcf.open(archaic_file);
  modern_vcf.open(modern_file);

  // write header
  output << "chrom\tpos\tref\talt";
  // build files, this will print the headers as well (archaic first)
  VCF_File archaic(&archaic_vcf, output);
  VCF_File modern(&modern_vcf, output);
  // terminate header
  output << "\n";

  bool recheck = false;
  // for each line in modern
  while (modern.update()) {
    // short update with recheck for case when archaic position is > modern
    while (recheck || archaic.update(true)) {
      recheck = false;

      // less than position, copy archaic and use modern blank line
      if (archaic.getPosition() < modern.getPosition()) {
        // skip lines with no informative archaic GT
        bool nonzero = false;
        auto gen = archaic.getGenotypes();
        for (int i = 0; i < 2 * archaic.getCount(); i += 2)
          if (gen[i] != '0') {
            nonzero = true;
            break;
          }

        // found at least one informative archaic site
        if (nonzero) {
          // output archaic information and blank modern
          output << archaic.getChromosome() << '\t' << archaic.getPosition()
                 << '\t' << archaic.getReference() << '\t'
                 << archaic.getAlternative() << '\t';
          output << archaic.getGenotypes();
          output << modern.getBlank();
          output << "\n";
        }
        // equal, have to check other conditions to write
      } else if (archaic.getPosition() == modern.getPosition()) {
        // check alleles are matching
        if (modern.isValid() &&
            archaic.getReference() == modern.getReference() &&
            (archaic.getAlternative() == '.' ||
             archaic.getAlternative() == modern.getAlternative())) {
          output << archaic.getChromosome() << '\t' << archaic.getPosition()
                 << '\t' << archaic.getReference() << '\t'
                 << modern.getAlternative() << '\t';
          output << archaic.getGenotypes();
          output << modern.getGenotypes();
          output << "\n";
        }
        // advance both to match legacy version
        break;
        // greater than, advance modern but keep archaic where it is
      } else {
        recheck = true;
        break;
      }
    }
  }

  if (of.is_open()) of.close();
  archaic_vcf.close();
  modern_vcf.close();

  return 0;
}
