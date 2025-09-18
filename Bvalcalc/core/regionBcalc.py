from Bvalcalc.core.chromBcalc import chromBcalc
from Bvalcalc.utils.load_bed_gff import load_bed_gff
from Bvalcalc.utils.load_Bmap import load_Bmap
from Bvalcalc.utils.header_utils import create_header_info_from_args, generate_header, write_headers_to_file
from Bvalcalc.core.calculateB import calculateB_unlinked
import numpy as np

def regionBcalc(args, region):    

    allblockstart, allblockend, allblockchrom  = load_bed_gff(args.bedgff) # Read BED/GFF, return start and end of conserved elements

    import Bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    dfe_helper.CONSTANT_DFE = args.constant_dfe # Update DFE if --constant_dfe

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

    # Create and write header once at the start
    if args.out is not None and not getattr(args, 'no_header', False):
        header_info = create_header_info_from_args(args, "B-map")
        header_lines = generate_header(header_info)
        write_headers_to_file(args.out, header_lines, mode='w')
    
    output_data, block_ranges, rec_rate_per_chunk_in_region, chunk_size = chromBcalc(args, blockstart, blockend, chromosome, unlinked_B, prior_pos, prior_b, calc_start, calc_end, caller="regionBcalc", write_header=False)

    return  output_data, block_ranges, rec_rate_per_chunk_in_region, chunk_size

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