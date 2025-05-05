from core.chromBcalc import chromBcalc
from core.utils.bedgffHandler import bedgffHandler
from core.utils.BmapHandler import BmapHandler
from core.calculateB import calculateB_unlinked
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
    blockstart, blockend = allblockstart[mask], allblockend[mask]
    chromosome = calc_chrom
    unlinked_blockstart, unlinked_blockend = allblockstart[~mask], allblockend[~mask]
    unlinked_L = np.sum(unlinked_blockend-unlinked_blockstart)
    unlinked_B = calculateB_unlinked(unlinked_L)

    if args.prior_Bmap is not None:
        prior_chromosomes, prior_positions, prior_b_values = BmapHandler(file_path = args.prior_Bmap)
        if not args.quiet: print(f"Using prior B values from {args.prior_Bmap}")
        prior_mask = (prior_chromosomes == chromosome)
        prior_pos = prior_positions[prior_mask]
        prior_b = prior_b_values[prior_mask]
    else:
        prior_pos, prior_b = None, None

    if args.out is not None: # Overwrite existing file with header
        with open(args.out, 'w') as out_f:
            out_f.write("Chromosome,Position,Conserved,B\n")
    output_data, block_ranges = chromBcalc(args, blockstart, blockend, chromosome, unlinked_B, prior_pos, prior_b, calc_start, calc_end, caller="regionBcalc")

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