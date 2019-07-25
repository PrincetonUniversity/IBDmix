# IBDmix
> Detecting archaic sequence segregation in modern human genomes

## Background
Admixture has played a prominent role in shaping patterns of human genomic
variation, including gene flow with now extinct hominins like Neanderthals and
Denisovans. We describe a novel probabilistic method called IBDmix to identify
introgressed hominin sequences, which unlike existing approaches, does not use
a modern reference population. 

## Usage

IBDmix consists of 3 steps:
1. Generate custom genotype file
2. Detect introgressed regions
3. Sort and filter results

### Compilation
#### Generate Genotype
`generate_gt` has a single source file and can be compiled with
```
g++ -std=c++11 IBDmix/generate_gt.cpp -o generate_gt
```

#### IBDmix
The `ibdmix` function requires several source files and is compiled with
```
g++ -std=c++11 IBDmix/IBDmix.cpp IBDmix/Genotype_Reader.cpp \
    IBDmix/IBD_Collection.cpp IBDmix/IBD_Segment.cpp IBDmix/IBD_Stack.cpp \
    -o ibdmix
```

#### Summary
The summary.sh script is a wrapper around awk and sort and
requires no compilation.

### Usage Details
#### Generate Genotype
`generate_gt` has the following options:
- __-h, --help__:    Print the help information and exit.
- __-a, --archaic__: The archaic vcf file.  May contain multiple
                 archaic samples but only one will be utilized for IBDmix.
                 Must be uncompressed text.
- __-m, --modern__:  The modern vcf file. Must be uncompressed text.
- __-o, --output__:  The merged genotype file output.  Written as uncompressed
                 text.
All files must be specified to run.

#### IBDmix
`ibdmix` takes the following options:
- __-h, --help__:       Print the help information and exit.
- __-g, --genotype__:   The input genotype file produced by `generate_gt`.
                        Required.
- __-o, --output__:     The output file.  Format is tab-delimited text with
                        columns for individual ID, chromosome, start, end, 
                        and LOD score.  The regions are half open with [start,
                        end), where the end position is the next position in
                        the input genotype file. Required.
- __-s, --sample__:     File specifying which samples to consider.  One
                        individual per line, must match header in genotype
                        file exactly.  If not specified all samples will be
                        utilized (all columns *except* archaic).
- __-n, --archaic__:    Name of archaic individual.  Must match column name.
                        If not specified, taken as the first column of the
                        genotype file (default output order for `generate_gt`).
- __-r, --mask__:       Bed file of masked regions.  If not specified all sites
                        are considered.
- __-d, --LOD-threahold__:       
                        Threshold value of log(odds) for emitting regions.
                        Default: 3.0
- __-m, --minor-allele-count-threshold__:       
                        Threshold count for filtering minor alleles.  For
                        each site, the number of occurances in the sample 
                        population is tabulated.  If fewer than this count
                        are found, the LOD score for the site is taken as 0.
                        Default: 1
- __-a, --archaic-error__:       
                        Allele error rate for archaic DNA.
                        Default: 0.01
- __-e, --modern-error-max__:       
                        Maximum allele error rate for modern humans.
                        Default: 0.0025
- __-c, --modern-error-proportion__:       
                        Ratio between allele error rate and minor allele
                        frequency for modern humans.  The error rate for a
                        given site is taken as the minimum of the `-e` option
                        and the product of this value with the allele
                        frequency.
                        Defualt: 2

#### Summary
Once a run of `ibdmix` completes, it is informative to filter the results
on a range of LOD values and length cutoffs.  It is faster to perform this
operation on the `ibdmix` output than to rerun with different options.

`summary.sh` takes up to four options in order:
- length cutoff:        The minimum length to emit a region. Use 0 for all.
- LOD cutoff:           The minimum LOD score to emit a region. Use 0 for all.
- input:                Optional, the uncompressed input file.
- output:               Optional, the uncompressed output file. Contains only
                        the regions passing both filters, sorted by ID then
                        start position of the region.

If neither input nor output is specified, standard input and output are 
utilized.  To specify only one value, use `-` to indicate standard 
input/output.  For example, to filter length greater than 1000 and LOD greater
than 5 from ibd\_output.txt and generate compressed output.gz
```
summary.sh 1000 5 ibd_output.txt - | gzip > output.gz
```

### Snakemake Workflow
The above workflow is automated through
[snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
In the provided snakefiles/Snakefile, all steps from compilation to summary
generation are encoded.  By default, summary files will be generated but
the workflow can be customized by setting options in snakefiles/config.yaml.
All files are compressed with gzip.

The easiest way to get started is to set your paths in the config file and
run `snakemake` in the snakefiles directory.

#### Important Configuration Options
> A note on wildcards.
>
> Wildcards are specified in snakemake by enclosing them in curly braces. When
changing options in the config file, it is important that all wildcards are
kept and spelled correctly.  For example, say you want to change the summary
output, which is currently `altai_1kg_{population}_{chrom}.gz`.  As long
as {chrom} and {population} are present, the rest can change.  If you want
to place each population in a separate directory, this could be changed to
`{population}/altai_1kg_chr_{chrom}.txt.gz`.  The output files would look like
`ASN/altai_1kg_chr_20.txt.gz`.  However, `{pop}/altai_1kg_{chr}.gz` will throw
and error because snakemake expects to find specific wildcards.

- compiler and options: These are utilized to make the executables from cpp
source code.  The compiler requires C++ 11 standard.
- exe\_root: Where to place compiled executables.
- source\_root: The path to the IBDmix directory.
- archaic\_vcf, modern\_vcf: Required inputs
- genotype\_file, ibd\_output: Temporary intermediate files
- ibd\_summary: Final output.  If no summary values are supplied, ibd\_output
is the final target instead.
- sample\_file: The input sample files for `ibdmix`. If removed, all samples
will be used and 'ALL' will be the population name.
- mask\_file: The bed file to use for each chromosome.  If removed, no
masking will be performed during IBDmix.
- IBDmix options: All parameters for performing ibdmix.  The genotype, output,
mask, and sample are handled automatically.  Another options can be supplied
here (including archaic sample name).
- IBDmix summary\_lod: List of all LOD values for filter.  If sorted output
is desired, a value of 0 will keep all regions. Remove this entry to just
produce the raw IBD files.
- IBDmix summary\_length: List of all length values for filter.  If sorted output
is desired, a value of 0 will keep all regions.  All permutations of length and
lod are produced.


#### Useful Snakemake Options
- --notemp: Without this option, genotype files and raw IBD calls will be
automatically deleted when they are no longer needed by the workflow.  If you
expect to perform multiple executions utilize this flag to keep those files
around.
- --delete-temp-output: Once you are done with temp files this will find and
delete all specified by the current config.yaml
- --cores: Specify the number of threads available to snakemake.  With more
cores, more jobs can be run simultaneously.
- --dryrun: See what rules will run without running any of them.

## License
GNU GPLv3
