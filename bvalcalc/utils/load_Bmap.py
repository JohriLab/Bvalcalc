import sys
import csv
import numpy as np

def load_Bmap(file_path):
    chromosomes = []
    positions   = []
    b_values    = []

    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            # skip comments or stray headers
            if row[0].startswith('#') or row[0] == 'Chromosome':
                continue
            if len(row) < 3:
                # malformed line, skip
                continue

            chrom = row[0]
            pos   = int(row[1])
            # row[2] is "Conserved" but we don't need it here
            b     = float(row[2])

            chromosomes.append(chrom)
            positions.append(pos)
            b_values.append(b)
        chromosomes = np.array(chromosomes, dtype='<U20')
        positions = np.array(positions,   dtype=np.int64)
        b_values = np.array(b_values,    dtype=np.float64)

    return chromosomes, positions, b_values