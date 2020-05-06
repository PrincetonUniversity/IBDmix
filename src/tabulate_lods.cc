#include <string.h>
#include <stdio.h>

#include <iostream>
#include <fstream>

#include <CLI/CLI.hpp>

#include "IBDmix/Genotype_Reader.h"

int main(int argc, char *argv[]) {
    CLI::App app{"Calculate population specific LOD scores for all sites"};

    std::string genotype_file;
    app.add_option("-g,--genotype", genotype_file, "The genotype file")
        ->check(CLI::ExistingFile)
        ->required();

    std::string outfile = "-";
    app.add_option("-o,--output", outfile, "The output file location");

    std::string sample_file = "";
    app.add_option("-s,--sample", sample_file, "File containing samples to "
            "select from genotype.  Default to all samples in genotype.")
        ->check(CLI::ExistingFile);

    std::string archaic = "";
    app.add_option("-n,--archaic", archaic, "Name of archaic sample, default"
            " to first sample in genotype file");

    std::string mask_file = "";
    app.add_option("-r,--mask", mask_file, "Mask of regions to 'remove'. "
            "Regions in bed file have LOD set to 0")
        ->check(CLI::ExistingFile);

    int ma_threshold = 1;
    double archaic_error = 0.01,
           modern_error_max = 0.0025,
           modern_error_prop = 2;
    bool include_ninfs = false;
    bool include_zeros = false;
    app.add_option("-m,--minor-allele-count-threshold", ma_threshold,
            "Threshold count for filtering minor alleles");
    app.add_option("-a,--archaic-error", archaic_error,
            "Allele error rate for archaic DNA");
    app.add_option("-e,--modern-error-max", modern_error_max,
            "Maximum allele error rate for modern samples");
    app.add_option("-c,--modern-error-proportion", modern_error_prop,
            "Ratio between allele error rate and minor allele frequency");
    app.add_flag("--include-ninfs", include_ninfs,
            "Include sites with -ninf as LOD scores.  Will translate to "
            "-100 as the LOD score");
    app.add_flag("--include-zeros", include_zeros,
            "Include sites where all LODs are 0. Commonly occurs for sites in "
            "masked regions.");

    CLI11_PARSE(app, argc, argv);

    std::ifstream genotype;
    genotype.open(genotype_file);

    std::ifstream sample;
    if (sample_file != "")
        sample.open(sample_file);
    std::ifstream mask;
    if (mask_file != "")
        mask.open(mask_file);

    std::ofstream of;
    std::streambuf * buf;
    if (outfile == "-") {
        buf = std::cout.rdbuf();
    } else {
        of.open(outfile);
        buf = of.rdbuf();
    }
    std::ostream output(buf);

    // write header
    output << "chrom\tpos\tref\talt\tarchaic\tfreq_b\t0\t1\t2\n";

    Genotype_Reader reader(
            &genotype,
            &mask,
            archaic_error,
            modern_error_max,
            modern_error_prop,
            1e-200,  // minesp
            ma_threshold);

    int num_samples = reader.initialize(sample, archaic);
    if (sample.is_open())
        sample.close();

    while (reader.update()) {
        // all sites are 0
        if ( !include_zeros && reader.lod_cache[0] == 0 &&
                reader.lod_cache[1] == 0 &&
                reader.lod_cache[2] == 0)
            continue;

        // any sites are inf, replace with -100
        if (isinf(reader.lod_cache[0])) {
            reader.lod_cache[0] = -100;
            if (!include_ninfs) continue;
        }

        if (isinf(reader.lod_cache[1])) {
            reader.lod_cache[1] = -100;
            if (!include_ninfs) continue;
        }

        if (isinf(reader.lod_cache[2])) {
            reader.lod_cache[2] = -100;
            if (!include_ninfs) continue;
        }

        output << reader.chromosome << '\t'
            << reader.position << '\t'
            << reader.ref << '\t'
            << reader.alt << '\t'
            << reader.archaic << '\t'
            << reader.allele_frequency << '\t'
            << reader.lod_cache[0] << '\t'
            << reader.lod_cache[1] << '\t'
            << reader.lod_cache[2] << '\n';
    }

    genotype.close();
    if (mask.is_open())
        mask.close();
    if (of.is_open())
        of.close();
}
