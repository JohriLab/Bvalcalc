import csv
import os
import numpy as np
import sys

def load_bed_gff(file_path):
    """
    Read a BED, CSV, GFF, GFF3, or GTF file and return arrays of block starts, block ends, and chromosomes.

    - BED: tab-delimited, cols: chrom, start, end
    - CSV: comma-delimited, cols: chrom, start, end
    - GFF/GFF3/GTF: tab-delimited, cols: chrom, source, feature, start, end, ...

    Any lines beginning with '#' or with fewer than the required columns are skipped.
    Duplicate (start, end) pairs are removed.
    """
    # Determine extension and parsing rules
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    if ext == '.bed':
        delim = '\t'
        chrom_idx, start_idx, end_idx = 0, 1, 2
        min_cols = 3
    elif ext in ('.gff', '.gff3', '.gtf'):
        delim = '\t'
        chrom_idx, start_idx, end_idx = 0, 3, 4
        min_cols = 5
    else:
        # default to CSV
        delim = ','
        chrom_idx, start_idx, end_idx = 0, 1, 2
        min_cols = 3

    blockstart = []
    blockend = []
    chromosomes = []
    seen_blocks = set()

    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=delim)
        for row in reader:
            # skip comments
            if row and row[0].startswith('#'):
                continue
            # skip too-short rows
            if len(row) < min_cols:
                continue

            chrom = row[chrom_idx]
            try:
                start = int(row[start_idx])
                end = int(row[end_idx])
            except ValueError:
                # skip rows with non-integer coordinates
                continue

            # ensure start <= end
            if end < start:
                start, end = end, start
            # if start == end: # Remove 0-length exons
            #     continue
            if start == end: # skip zero-length elements
                continue

            # deduplicate
            if (start, end) not in seen_blocks:
                seen_blocks.add((start, end))
                chromosomes.append(chrom)
                blockstart.append(start)
                blockend.append(end)

    return np.array(blockstart), np.array(blockend), np.array(chromosomes)
