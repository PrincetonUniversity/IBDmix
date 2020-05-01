# translate physical distances in ibd_bed file to genetic distance
import gzip
import sys


class gen_mapper:
    def __init__(self, reader):
        self.reader = reader
        self.bp_positions = []
        self.cm_positions = []
        self.finished = False
        gen_map.readline()  # remove header
        self.read_interval()
        self.read_interval()

    def read_interval(self):
        line = self.reader.readline()
        # if reader is empty
        if line == "":
            self.finished = True
            return

        toks = line.split()
        bp = float(toks[1])
        cm = float(toks[3])
        self.bp_positions.append(bp)
        self.cm_positions.append(cm)

    def predict(self, position):
        if position < self.bp_positions[0]:
            return self.cm_positions[0]

        while position > self.bp_positions[-1]:
            if not self.finished:
                # read new line
                self.read_interval()
            else:
                return self.cm_positions[-1]

        # find location with neighboring positions, will be closest to the
        # back so will iterate in reverse
        for i, el in enumerate(reversed(self.bp_positions)):
            if i == 0:
                continue
            ind = len(self.bp_positions) - i - 1
            el2 = self.bp_positions[ind + 1]
            if el <= position < el2:
                break

        return ((position - self.bp_positions[ind]) *
                ((self.cm_positions[ind] - self.cm_positions[ind + 1]) /
                 (self.bp_positions[ind] - self.bp_positions[ind + 1]))
                + self.cm_positions[ind])

if 'snakemake' in locals() or 'snakemake' in globals():
    gen = snakemake.input[0]
    ibd = snakemake.input[1]
    output = snakemake.output[0]
else:
    gen = sys.argv[1]
    ibd = sys.argv[2]
    output = sys.argv[3]

with open(gen, 'r') as gen_map, \
        gzip.open(ibd, 'rt') as ibd, \
        gzip.open(output, 'wt') as output:
    mapper = gen_mapper(gen_map)
    output.write(ibd.readline().strip() + '\tcm_dist\n')  # copy header
    for line in ibd:
        toks = line.split()
        start = int(toks[2])
        end = int(toks[3])
        cm_dist = mapper.predict(end) - mapper.predict(start)
        output.write(line.strip() + f'\t{cm_dist:0.7f}\n')
