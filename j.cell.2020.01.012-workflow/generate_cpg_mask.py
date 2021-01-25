import gzip
import re
from collections import defaultdict


def main(infiles, outfile):
    snps = {species: simple_snps(infiles)
            for species, infiles in infiles.items()
            if species != 'reference' and species != 'modern'}

    snps['modern'] = modern_snps(infiles['modern'])

    skipping = False

    for line in open(infiles['reference'], 'r'):
        if line.startswith('>'):
            pos = 1
            this_base_C = False
            this_base_G = False
            last_base_C = False

            chrom = re.sub('^>', '', line.split()[0])
            if chrom.isnumeric():
                skipping = False
                for snp in snps.values():
                    snp.next_chrom(chrom)
            else:
                skipping = True
            continue

        if not skipping:
            line = line.strip().upper()

            for i, base in enumerate(line):
                this_base_C = base == 'C'
                this_base_G = base == 'G'
                for snp in snps.values():
                    C, G = snp.is_cpg(pos+i, base)
                    this_base_C |= C
                    this_base_G |= G

                if last_base_C and this_base_G:
                    outfile.write(f'{chrom}\t{pos-2+i}\t{pos+i}\n')

                last_base_C = this_base_C

            pos += len(line)


class simple_snps:
    def __init__(self, infiles):
        self.infiles = {
            re.split('[-_.]', file)[-2]: file
            for file in infiles}
        self.infile = None
        self.position = None
        self.ref = None
        self.variant = None
        self.chrom = None

    def next_chrom(self, chrom):
        if self.infile is not None:
            self.infile.close()

        self.infile = gzip.open(self.infiles[chrom], 'rt')
        self.update()

    def update(self):
        line = self.infile.readline().strip()
        if line == '':
            self.position = None
        else:
            self.chrom, self.position, self.ref, self.variant = line.split()
            self.position = int(self.position)

    def is_cpg(self, position, base):
        while self.position is not None and position > self.position:
            self.update()
        if position == self.position:
            assert (self.ref == base or self.ref in ('N', 'M')
                    or base in ('N', 'M')), (f'Expected {position} '
                                             f'{base} to match {self.ref}')
            return self.variant == 'C', self.variant == 'G'
        else:  # position < self.position
            return False, False


class modern_snps:
    # modern snps have all chroms in one file, read all then access as needed
    def __init__(self, infile):
        self.data = defaultdict(dict)
        with gzip.open(infile, 'rt') as infile:
            for line in infile:
                chrom, pos, ref, alt = line.split()
                pos = int(pos)
                self.data[chrom][pos] = (ref, alt)
        self.current_chrom = None

    def next_chrom(self, chrom):
        self.current_chrom = self.data[chrom]

    def is_cpg(self, position, base):
        if position in self.current_chrom:
            ref, alt = self.current_chrom[position]
            assert (ref == base or ref in ('N', 'M')
                    or base in ('N', 'M')), (f'Expected {position} '
                                             f'{base} to match {ref}')
            return alt == 'C', alt == 'G'
        else:  # position not found
            return False, False


if __name__ == '__main__':
    infiles = snakemake.input

    outfile = snakemake.output[0]
    with open(outfile, 'w') as output:
        main(infiles, output)
