import csv
import numpy as np

def bedgffHandler(file_path):
    blockstart = []
    blockend = []
    chromosomes = []  # List to hold the chromosome name for each block
    seen_blocks = set()

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row and row[0].startswith("#"):  # Skip header or comment lines
                continue
            if len(row) >= 3:
                chrom = str(row[0])
                start, end = int(row[1]), int(row[2])

                if end < start:
                    start, end = end, start  # Flip if needed

                if (start, end) not in seen_blocks:
                    seen_blocks.add((start, end))
                    chromosomes.append(chrom)
                    blockstart.append(start)
                    blockend.append(end)

    blockstart = np.array(blockstart)
    blockend = np.array(blockend)
    chromosomes = np.array(chromosomes)

    return blockstart, blockend, chromosomes
