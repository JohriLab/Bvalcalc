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
        print(blockstart, blockend, chromosome)
        output_data, block_ranges = chromBcalc(args, blockstart, blockend, chromosome)
        print("First chromosomes output done, need to figure out what to do with it now...")
        sys.exit()

    # GET UNLINKED B FOR EACH CHROMOSOME
    # FOR EACH UNIQUE CHROM
    # # DO: 
    
    # # # NO INPUT COLLECTED!!
    # # Pass up plot data as needed

    return output_data, block_ranges