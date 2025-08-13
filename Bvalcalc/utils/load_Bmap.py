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

            # Strip trailing single quote if present
            b_str = row[2].strip()
            if b_str.endswith("'"):
                b_str = b_str[:-1]
            b = float(b_str)

            chromosomes.append(chrom)
            positions.append(pos)
            b_values.append(b)

    chromosomes = np.array(chromosomes, dtype='<U20')
    positions   = np.array(positions,   dtype=np.int64)
    b_values    = np.array(b_values,    dtype=np.float64)

    return chromosomes, positions, b_values
