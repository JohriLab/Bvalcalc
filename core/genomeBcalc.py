from core.chromBcalc import chromBcalc
from core.utils.bedgffHandler import bedgffHandler
import numpy as np

def genomeBcalc(args):    

    blockstart, blockend, blockchrom = bedgffHandler(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements

    unique_chroms = np.unique(blockchrom) # Move BED/GFF handler here
    # GET LIST OF EACH UNIQUE CHROMOSOME
    # GET UNLINKED B FOR EACH CHROMOSOME
    # FOR EACH UNIQUE CHROM
    # # DO: 
    output_data, block_ranges = chromBcalc(args, blockstart, blockend, blockchrom)
    # # # NO INPUT COLLECTED!!
    # # Pass up plot data as needed

    return output_data, block_ranges