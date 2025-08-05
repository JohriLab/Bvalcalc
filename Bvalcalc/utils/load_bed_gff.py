import csv
import os
import numpy as np
import re

def load_bed_gff(file_path):
    """
    Read a BED, CSV, GFF, GFF3, or GTF file and return arrays of block starts, block ends, and chromosomes.

    - BED: tab-delimited, cols: chrom, start, end
    - CSV: comma-delimited, cols: chrom, start, end
    - GFF/GFF3/GTF: tab-delimited, cols: chrom, source, feature, start, end, ...

    Any lines beginning with '#' or with fewer than the required columns are skipped.
    For duplicate start coordinates, only the longest block is kept.
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
        min_cols = 9
    else:
        # default to CSV
        delim = ','
        chrom_idx, start_idx, end_idx = 0, 1, 2
        min_cols = 3

    # mapping from (chrom, start) -> (end, length)
    longest_blocks = {}
    skip_count = 0

    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=delim)
        for row in reader:
            # skip comments
            if row and row[0].startswith('#'):
                continue
            # skip too-short rows
            if len(row) < min_cols:
                skip_count += 1
                continue

            chrom = row[chrom_idx]
            try:
                start = int(row[start_idx])
                end = int(row[end_idx])
            except ValueError:
                continue

            if end < start:
                start, end = end, start
            if start == end:
                continue

            key = (chrom, start)
            length = end - start

            if key not in longest_blocks or length > (longest_blocks[key][0] - start):
                longest_blocks[key] = (end, chrom)

    if skip_count:
        raise ValueError(f"ERROR: {skip_count} malformed lines (fewer than {min_cols} columns) in {file_path}. Consider re-formatting to .BED")

    # Unpack results
    # Convert dict to list of (chrom, start, end)
    raw_blocks = [(chrom, start, end)
                  for (chrom, start), (end, _) in longest_blocks.items()]

    # Sort by chrom, start, then descending end
    raw_blocks.sort(key=lambda x: (x[0], x[1], -x[2]))

    # Remove overlapping entries (keep longest non-overlapping)
    filtered = []
    current_end = {}

    for chrom, start, end in raw_blocks:
        if chrom in current_end and start < current_end[chrom]:
            continue
        filtered.append((chrom, start, end))
        current_end[chrom] = end

    # Unpack and return
    if filtered:
        chromosomes, blockstart, blockend = zip(*filtered)
        return np.array(blockstart), np.array(blockend), np.array(chromosomes)
    else:
        return np.array([]), np.array([]), np.array([])
