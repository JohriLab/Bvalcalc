from bvalcalc.utils.load_vcf import load_vcf
from bvalcalc.utils.load_Bmap import load_Bmap
import numpy as np


def vcfBmap(args, vcf_path):    

    vcf_chroms, vcf_pos = load_vcf(vcf_path)
    print(vcf_chroms, vcf_pos)

    bmap_chroms, bmap_pos, b_values = load_Bmap(file_path = args.Bmap)
    print(bmap_chroms, bmap_pos, b_values)

    max_bmap_pos = bmap_pos[-1]
    above_count  = np.sum(vcf_pos > max_bmap_pos)
    if above_count > 0:
        print(f"WARNING: {above_count} VCF positions are above the max B-map start position ({max_bmap_pos}) consider calculating an extended B-map")

    vcf_unique  = np.unique(vcf_chroms)
    bmap_unique = np.unique(bmap_chroms)
    missing_in_bmap = set(vcf_unique) - set(bmap_unique)
    missing_in_vcf  = set(bmap_unique) - set(vcf_unique)
    if missing_in_bmap:
        print(f"WARNING: The following chromosomes are in the VCF but not in the B-map: {missing_in_bmap}")
    if missing_in_vcf:
        print(f"WARNING: The following chromosomes are in the B-map but not in the VCF: {missing_in_vcf}")


    idx = np.searchsorted(bmap_pos, vcf_pos, side='right') - 1
    vcf_b_values = b_values[idx]



    # vcf_bvals = np.full(vcf_pos.shape, np.nan, dtype=np.float64)

    #     for chrom in bmap_chroms:
    #     mask_b = (bmap_chroms_raw == chrom)
    #     starts = bmap_pos[mask_b]
    #     vals   = bmap_vals[mask_b]

    print(idx, vcf_b_values)

    print("gottem in vcfBmap")

# 1. Return histogram of B for VCF positions
###     Also give B for each VCF position in file
# 2. Get position of sites with B above or below X, or top N positions of bottom N positions
###     Also recode VCF with those positions if specified

