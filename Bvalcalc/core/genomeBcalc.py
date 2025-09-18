from Bvalcalc.core.chromBcalc import chromBcalc
from Bvalcalc.utils.load_bed_gff import load_bed_gff
from Bvalcalc.utils.load_chr_sizes import load_chr_sizes
from Bvalcalc.utils.load_Bmap import load_Bmap
from Bvalcalc.utils.write_chrom_B_to_file import write_chrom_B_to_file
from Bvalcalc.utils.header_utils import create_header_info_from_args, generate_header, write_headers_to_file
from Bvalcalc.core.calculateB import calculateB_unlinked
import numpy as np
import re

def sort_chromosomes_naturally(chromosomes):
    """
    Sort chromosomes in natural order: chr1, chr2, chr3, ..., chr9, chr10, chr11, chrX, chrY, chrM
    """
    def natural_sort_key(chrom):
        # Extract the chromosome name part (remove 'chr' prefix if present)
        chrom_str = str(chrom)
        if chrom_str.startswith('chr'):
            chrom_str = chrom_str[3:]
        
        # Extract numeric part for natural sorting
        match = re.match(r'^(\d+)', chrom_str)
        if match:
            return int(match.group(1))
        
        # For non-numeric chromosomes, sort alphabetically
        return chrom_str
    
    return sorted(chromosomes, key=natural_sort_key)

def genomeBcalc(args):    

    allblockstart, allblockend, allblockchrom = load_bed_gff(args.bedgff) # Read BED/GFF, return start and end of conserved elements

    unique_chromosomes = np.unique(allblockchrom) # Move BED/GFF handler here
    unique_chromosomes = sort_chromosomes_naturally(unique_chromosomes)
    if args.chr_sizes is not None: 
        chr_sizes = load_chr_sizes(args.chr_sizes)  # <-- Path to the sizes CSV file

    import Bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    dfe_helper.CONSTANT_DFE = args.constant_dfe # Update DFE if --constant_dfe

    print("Chromosomes loaded:", unique_chromosomes) ## Now, loop over each chromosome and save B output

    if args.prior_Bmap is not None:
        prior_chromosomes, prior_positions, prior_b_values = load_Bmap(file_path = args.prior_Bmap)
        if not args.quiet: print(f"Using prior B values from {args.prior_Bmap}")

    # Create and write header once at the start
    if args.out is not None and not getattr(args, 'no_header', False):
        header_info = create_header_info_from_args(args, "B-map")
        header_lines = generate_header(header_info)
        write_headers_to_file(args.out, header_lines, mode='w')
    
    for i in np.arange(0,len(unique_chromosomes)):

        mask = allblockchrom == unique_chromosomes[i]
        blockstart, blockend = allblockstart[mask], allblockend[mask]
        chromosome = unique_chromosomes[i]
        unlinked_blockstart, unlinked_blockend = allblockstart[~mask], allblockend[~mask]
        unlinked_L = np.sum(unlinked_blockend-unlinked_blockstart)
        unlinked_B = calculateB_unlinked(unlinked_L)

        if args.prior_Bmap is not None: 
            prior_mask = (prior_chromosomes == chromosome)
            prior_pos = prior_positions[prior_mask]
            prior_b = prior_b_values[prior_mask]
        else:
            prior_pos, prior_b = None, None

        if args.chr_sizes is not None: chr_size = chr_sizes.get(chromosome)
        else: chr_size = None
        chromBcalc(args, blockstart, blockend, chromosome, unlinked_B, prior_pos, prior_b, calc_start=None, calc_end=None, chr_size=chr_size, caller="genomeBcalc", write_header=False)

    return