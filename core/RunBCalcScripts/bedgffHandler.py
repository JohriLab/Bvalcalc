import csv
import numpy as np

def bedgffHandler(file_path):

        ## 1. Read in BED input
    blockstart = []
    blockend = []
    seen_blocks = set()

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
                        if row and row[0].startswith("#"): # Skip header rows
                            continue
                        if len(row) >= 3:
                            start, end = int(row[1]), int(row[2])
                            if end < start:
                                start, end = end, start # Flip reverse orientation genes
                            if (start, end) not in seen_blocks:  # Check if the pair is unique
                                seen_blocks.add((start, end))  # Mark this pair as seen to avoid duplicates
                                blockstart.append(start)
                                blockend.append(end)
    blockstart = np.array(blockstart) # Position where blocks (selected regions) start
    blockend = np.array(blockend) # Position where blocks (selected regions) end

    return blockstart, blockend


