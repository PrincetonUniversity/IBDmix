import dask.dataframe as dd
from typing import TextIO
import os
import gzip


class File_Merger():
    '''
    Object to merge vcf and genotype files to a single genotype file
    '''
    def __init__(self, archaic_vcf: str, modern_vcf: str):
        self.archaic_vcf = archaic_vcf
        self.modern_vcf = modern_vcf

    def run(self, output: str):
        if os.path.splitext(self.archaic_vcf)[1] == '.gz':
            archaic_reader = gzip.open(self.archaic_vcf, 'rt')
        else:
            archaic_reader = open(self.archaic_vcf, 'r')
        archaic = VCF_File(archaic_reader, self.archaic_vcf)
        archaic_reader.close()

        if os.path.splitext(self.modern_vcf)[1] == '.gz':
            modern_reader = gzip.open(self.modern_vcf, 'rt')
        else:
            modern_reader = open(self.modern_vcf, 'r')
        modern = VCF_File(modern_reader, self.modern_vcf)
        modern_reader.close()

        result = self.merge(archaic, modern)

        result.to_csv(output,
                      sep='\t')

    def merge(self, archaic, modern) -> dd.DataFrame:
        '''
        Merge two dataframe chunks, returning the intersection
        '''
        result = dd.merge(archaic.get_dataframe(skip_non_informative=True),
                          modern.get_dataframe(),
                          on='POS', how='left',
                          suffixes=('_arch', '_mod'))
        result = result.fillna(0)
        # TODO either filter 9's here or when reading below
        result = result.loc[
            ((result['CHROM_mod'] == 0) &
             (result[archaic.get_sample_header()].sum(axis=1) != 0)) |
            ((result['CHROM_mod'] != 0) &
             (result['REF_arch'] == result['REF_mod']) &
             ((result['ALT_arch'] == result['ALT_mod']) |
              (result['ALT_arch'] == '.')
              ))
        ]
        result['ALT_arch'] = result['ALT_arch'].where(
            result['CHROM_mod'] == 0,
            result['ALT_mod'])
        result = result.drop(['CHROM_mod', 'REF_mod', 'ALT_mod'],
                             axis='columns')
        result = result.rename(columns={
            'CHROM_arch': 'CHROM',
            'REF_arch': 'REF',
            'ALT_arch': 'ALT',
        })
        result = result.astype({name: 'int32'
                                for name in modern.get_sample_header()})
        return result


class VCF_File():
    '''
    Object managing the import of vcf files as pandas tables
    '''
    def __init__(self, reader: TextIO, filename: str):
        self.reader = reader
        # get header stuff
        for line in reader:
            if line[1] != '#':
                self.header = line[1:].split()
                break

        self.filename = filename

    def get_dataframe(self, skip_non_informative=False):
        '''
        Return formatted genotype table
        '''
        gt_values = {
            '00': 0,
            '01': 1,
            '10': 1,
            '11': 2,
            '..': 9
        }
        converters = {h: lambda entry: gt_values[entry[0] + entry[2]]
                      for h in self.get_sample_header()}

        result = dd.read_csv(self.filename,
                             sep='\t',
                             names=self.header,
                             usecols=([0, 1, 3, 4] +
                                      list(range(9, len(self.header)))),
                             comment='#',
                             converters=converters).set_index('POS',
                                                              sorted=True)
        result = result.loc[(result['REF'].str.len() == 1) &
                            (result['ALT'].str.len() == 1)]

        # remove sites where all entires are 9
        if skip_non_informative:
            result = result.loc[(result[self.get_sample_header()]
                                 != 9).any(axis=1)]

        return result

    def get_initial_header(self):
        return self.header[0:2] + self.header[3:5]

    def get_sample_header(self):
        return self.header[9:]
