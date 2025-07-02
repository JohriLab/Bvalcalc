from bvalcalc.core.chromBcalc import chromBcalc
from bvalcalc.utils.load_bed_gff import load_bed_gff
from bvalcalc.utils.load_chr_sizes import load_chr_sizes
from bvalcalc.utils.load_Bmap import load_Bmap
from bvalcalc.core.calculateB import calculateB_unlinked
import numpy as np

def genomeBcalc(args):    

    allblockstart, allblockend, allblockchrom = load_bed_gff(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements

    unique_chromosomes = np.unique(allblockchrom) # Move BED/GFF handler here
    if args.chr_sizes is not None: 
        chr_sizes = load_chr_sizes(args.chr_sizes)  # <-- Path to the sizes CSV file

    import bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe

    print("Chromosomes loaded:", unique_chromosomes) ## Now, loop over each chromosome and save B output

    if args.prior_Bmap is not None:
        prior_chromosomes, prior_positions, prior_b_values = load_Bmap(file_path = args.prior_Bmap)
        if not args.quiet: print(f"Using prior B values from {args.prior_Bmap}")

    if args.out is not None: # Overwrite existing file with header
        with open(args.out, 'w') as out_f:
            out_f.write("Chromosome,Position,B\n")

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
        chromBcalc(args, blockstart, blockend, chromosome, unlinked_B, prior_pos, prior_b, calc_start=None, calc_end=None, chr_size=chr_size, caller="genomeBcalc")

    return