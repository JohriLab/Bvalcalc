from core.chromBcalc import chromBcalc
import numpy as np

def genomeBcalc(args):    

    unique_chroms = np.unique() # Move BED/GFF handler here
    # GET LIST OF EACH UNIQUE CHROMOSOME
    # GET UNLINKED B FOR EACH CHROMOSOME
    # FOR EACH UNIQUE CHROM
    # # DO: 
    output_data, block_ranges = chromBcalc(args)
    # # # NO INPUT COLLECTED!!
    # # Pass up plot data as needed

    return output_data, block_ranges