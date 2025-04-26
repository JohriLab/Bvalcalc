from core.chromBcalc import chromBcalc
from core.utils.bedgffHandler import bedgffHandler
from core.utils.chrSizesHandler import load_chr_sizes
import numpy as np
import sys

def genomeBcalc(args):    

    allblockstart, allblockend, allblockchrom = bedgffHandler(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements

    unique_chromosomes = np.unique(allblockchrom) # Move BED/GFF handler here
    if args.chr_sizes is not None: 
        chr_sizes = load_chr_sizes(args.chr_sizes)  # <-- Path to the sizes CSV file

    import core.utils.dfeHelper as dfeHelper
    dfeHelper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe

    print("Chromosomes loaded:", unique_chromosomes) ## Now, loop over each chromosome and save B output

    if args.out is not None: # Overwrite existing file with header
        with open(args.out, 'w') as out_f:
            out_f.write("Chromosome,Position,Conserved,B\n")

    for i in np.arange(0,len(unique_chromosomes)):
        mask = allblockchrom == unique_chromosomes[i]
        blockstart = allblockstart[mask]
        blockend = allblockend[mask]
        chromosome = unique_chromosomes[i]
        if args.chr_sizes is not None: chr_size = chr_sizes.get(chromosome)
        else: chr_size = None
        chromBcalc(args, blockstart, blockend, chromosome, calc_start=None, calc_end=None, chr_size=chr_size, caller="genomeBcalc")

    return