from core.chromBcalc import chromBcalc
from core.utils.bedgffHandler import bedgffHandler
from core.utils.chrSizesHandler import load_chr_sizes
import numpy as np
import sys

def genomeBcalc(args):    

    allblockstart, allblockend, allblockchrom = bedgffHandler(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements

    unique_chromosomes = np.unique(allblockchrom) # Move BED/GFF handler here
    chr_sizes = load_chr_sizes(args.chr_sizes)  # <-- Path to the sizes CSV file

    print(unique_chromosomes) ## Now, loop over each chromosome and save B output
    for i in np.arange(0,len(unique_chromosomes)):
        mask = allblockchrom == unique_chromosomes[i]
        blockstart = allblockstart[mask]
        blockend = allblockend[mask]
        chromosome = unique_chromosomes[i]
        chr_size = chr_sizes.get(chromosome)
        chromBcalc(args, blockstart, blockend, chromosome, chr_size, caller="genomeBcalc")
    
    # # Don't need to capture output, simply save in chromBcalc or not. 
    # # Then, if you want to plot a region, must use --calc_region chr:start-end (instead of calc_start,calc_end)
    # # Pass up plot data as needed

    return