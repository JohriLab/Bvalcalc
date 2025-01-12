import csv
import numpy as np
from helperScripts.findFlankLen import findFlankLen 

def bedgffHandler(file_path, ):

        ## 1. Read in BED input
    blockstart = []
    blockend = []
    seen_blocks = set()

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 2:
                start, end = int(row[1]), int(row[2])
                if (start, end) not in seen_blocks:  # Check if the pair is unique
                    seen_blocks.add((start, end))  # Mark this pair as seen to avoid duplicates
                    blockstart.append(start)
                    blockend.append(end)
    blockstart = np.array(blockstart)
    blockend = np.array(blockend)
    lengths = blockend - blockstart #XX Need to extend for reverse orientation

# Calculate relevant flanking distances for each block (gene)
    flank_distances = np.zeros_like(lengths, dtype=np.int32)
    flank_blockstart = np.zeros_like(blockstart, dtype=np.int32)
    flank_blockend = np.zeros_like(blockend, dtype=np.int32)
    for i, length in enumerate(lengths):
        flank_distances[i] = findFlankLen(0.998, length)
        flank_blockstart[i] = blockstart[i] - flank_distances[i]
        flank_blockend[i] = blockend[i] + flank_distances[i]

    return blockstart, blockend, lengths, flank_blockstart, flank_blockend 


