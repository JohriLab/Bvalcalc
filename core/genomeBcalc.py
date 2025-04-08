from core.chromBcalc import chromBcalc
from core.utils.bedgffHandler import bedgffHandler
import numpy as np
import sys

def genomeBcalc(args):    

    allblockstart, allblockend, allblockchrom = bedgffHandler(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements

    unique_chromosomes = np.unique(allblockchrom) # Move BED/GFF handler here
    print(unique_chromosomes) ## Now, loop over each chromosome and save B output
    for i in np.arange(0,len(unique_chromosomes)):
        mask = allblockchrom == unique_chromosomes[i]
        blockstart = allblockstart[mask]
        blockend = allblockend[mask]
        chromosome = unique_chromosomes[i]
        output_data, block_ranges = chromBcalc(args, blockstart, blockend, chromosome)
    
    # # Don't need to capture output, simply save in chromBcalc or not. 
    # # Then, if you want to plot a region, must use --calc_region chr:start-end (instead of calc_start,calc_end)
    # # Pass up plot data as needed

    return output_data, block_ranges