import click
from merge_genotype import File_Merger


@click.group()
def IBDmix():
    pass


@IBDmix.command()
@click.option('-m', '--modern-vcf',
              help='modern vcf file')
@click.option('-a', '--archaic-vcf',
              help='archaic vcf file')
@click.option('-o', '--output',
              help='Output file, txt or gz')
def merge_genotype(modern_vcf, archaic_vcf, output):
    '''
    Merge archaic and modern vcf files into one genotype file
    '''
    merger = File_Merger(archaic_vcf, modern_vcf)
    merger.run(output)


if __name__ == '__main__':
    IBDmix()
