from bvalcalc.core.chromBcalc import chromBcalc
from bvalcalc.utils.load_bed_gff import load_bed_gff
from bvalcalc.utils.load_Bmap import load_Bmap
from bvalcalc.core.calculateB import calculateB_unlinked
import numpy as np
import sys

def regionBcalc(args, region):    

    allblockstart, allblockend, allblockchrom,  = load_bed_gff(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements

    import bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe

    calc_chrom, calc_start, calc_end = parse_region(region)

    mask = allblockchrom == calc_chrom
    blockstart, blockend = allblockstart[mask], allblockend[mask]
    chromosome = calc_chrom
    unlinked_blockstart, unlinked_blockend = allblockstart[~mask], allblockend[~mask]
    unlinked_L = np.sum(unlinked_blockend-unlinked_blockstart)
    unlinked_B = calculateB_unlinked(unlinked_L)

    if args.prior_Bmap is not None:
        prior_chromosomes, prior_positions, prior_b_values = load_Bmap(file_path = args.prior_Bmap)
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