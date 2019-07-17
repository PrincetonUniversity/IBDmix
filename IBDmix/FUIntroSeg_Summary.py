import pandas as pd
import sys
import getopt


def remove_overlapSegements(out, sub, Anum):
    '''
    If multiple archaic are present, find max overlapping region
    Add in the size of the IBD segment as a final column
    '''
    prev_row = sub.iloc[0, :]
    with open(out, 'wt') as f:
        f.write('ID\tchr\tstart\tend\t')
        for k in range(Anum):
            f.write('Archaic_%d\t' % (k))
        f.write('MaxLOD\tsize\n')

        for index, row in sub.iterrows():
            if index >= 0 and row.ID == prev_row.ID \
                    and row.chr == prev_row.chr and row.start < prev_row.end:
                if row.MaxLOD > prev_row.MaxLOD:
                    prev_row = row
            else:
                f.write('%s\t%d\t%d\t%d\t' % (prev_row.ID, prev_row.chr,
                                              prev_row.start, prev_row.end))
                for k in range(Anum):
                    f.write('%g\t' % prev_row['Archaic_'+str(k)])
                f.write('%g\t%d\n' % (prev_row.MaxLOD, prev_row['size']))
                prev_row = row
        f.write('%s\t%d\t%d\t%d\t' % (prev_row.ID, prev_row.chr,
                                      prev_row.start, prev_row.end))
        for k in range(Anum):
            f.write('%g\t' % prev_row['Archaic_'+str(k)])
        f.write('%g\t%d\n' % (prev_row.MaxLOD, prev_row['size']))


def print_help():
    '''
    Print help message to stdout
    '''
    print('python ./FUIntroSeg_Summary.py -f <filename> -i <#Archaic Homoids> '
          '-d <cutoff for lod score> -l <cutoff for segment size>')
    print('-h, --help:           print this help and exit')
    print('-f, --file:           input file')
    print('-i, --num:            number of archaic homoids')
    print('-d, --lod:            cutoff for lod score')
    print('-l, --len:            cutoff for segment size')


try:
    opts, args = getopt.getopt(sys.argv[1:], 'hf:i:d:l:',
                               ['help', 'file=', 'num=', 'LOD=', 'Length='])
except getopt.GetoptError:
    print_help()
    sys.exit(1)
for opt, arg in opts:
    if opt in ('-h', '--help'):
        print_help()
        sys.exit(0)
    elif opt in ('-f', '--file'):
        filename = arg
    elif opt in ('-i', '--num'):
        Anum = int(arg)
    elif opt in ('-d', '--lod'):
        cutoff_lod = float(arg)
    elif opt in ('-l', '--len'):
        cutoff_len = int(arg)

dat = pd.read_csv(filename, sep='\t', header=None,
                  names=['ID', 'chr', 'start', 'end'] +
                  ['Archaic_'+str(k) for k in range(Anum)] + ['MaxLOD'])
dat['size'] = dat.end - dat.start
# filter data based on lod cutoff and size
dat = dat.loc[(dat['MaxLOD'] > cutoff_lod) & (dat['size'] > cutoff_len), :]
# sort
dat = dat.sort_values(['ID', 'chr', 'start', 'end'])
remove_overlapSegements(
    f'ALL_D{cutoff_lod}_L{cutoff_len}_{filename}', dat, Anum)

dat = pd.read_csv(
    f'ALL_D{cutoff_lod}_L{cutoff_len}_{filename}', sep='\t')
# add in sel_N for Nth archaic, where sel_N == 1 if it was the maxLOD
for k in range(Anum):
    dat['sel_'+str(k)] = dat.apply(
        lambda row: 1 if row['Archaic_'+str(k)] == row['MaxLOD'] else 0,
        axis=1)
lst_cols = ['sel_'+str(k) for k in range(Anum)]
# count how many archaics match maxLOD
dat['count'] = dat[lst_cols].sum(axis=1)

# write sites found in a single archaic (archaic{k}) or ambiguous
# (with ambiguous including the single archaic)
for k in range(Anum):
    sub = dat.loc[(dat['count'] == 1) & (dat['sel_'+str(k)] == 1),
                  ['ID', 'chr', 'start', 'end', 'MaxLOD']]
    sub.to_csv(
        f'./Merged_Archaic{k}_D{cutoff_lod}_L{cutoff_len}_{filename}',
        sep='\t', index=False)
    sub = dat.loc[(dat['sel_'+str(k)] == 1),
                  ['ID', 'chr', 'start', 'end', 'MaxLOD']]
    sub.to_csv(
        f'./MergedAmbig_Archaic{k}_D{cutoff_lod}_L{cutoff_len}_{filename}',
        sep='\t', index=False)
