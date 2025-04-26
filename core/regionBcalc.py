from core.chromBcalc import chromBcalc
from core.utils.bedgffHandler import bedgffHandler
import numpy as np
import sys

def regionBcalc(args):    

    allblockstart, allblockend, allblockchrom,  = bedgffHandler(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements
    import core.utils.dfeHelper as dfeHelper
    dfeHelper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe

    calc_chrom, calc_start, calc_end = parse_region(args.calc_region)

    # # First, require calc_region
    # # Then, mask for the chromosome the calc_region is in.
    # # Then, pass that into chromBcalc() as calc_start and calc_end
    mask = allblockchrom == calc_chrom
    blockstart = allblockstart[mask]
    blockend = allblockend[mask]
    chromosome = calc_chrom
    if args.out is not None: # Overwrite existing file with header
        with open(args.out, 'w') as out_f:
            out_f.write("Chromosome,Position,Conserved,B\n")
    output_data, block_ranges = chromBcalc(args, blockstart, blockend, chromosome, calc_start, calc_end, caller="regionBcalc")

    return  output_data, block_ranges

def parse_region(region_str):
    try:
        chrom_part, pos_part = region_str.split(":")
        start_str, end_str = pos_part.split("-")
        chrom = chrom_part
        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))
        return chrom, start, end
    except ValueError:
        raise ValueError(f"Region format invalid: '{region_str}' should be like 'chr1:12345-67890'")