from bvalcalc.utils.load_vcf import load_vcf
from bvalcalc.utils.load_Bmap import load_Bmap
import numpy as np
import csv
import sys

def vcfBmap(args, vcf_path):    
    #Load info using utils
    vcf_chroms, vcf_pos = load_vcf(vcf_path)
    bmap_chroms, bmap_pos, b_values = load_Bmap(file_path=args.Bmap)

    max_bmap_pos = bmap_pos[-1]
    above_count = np.sum(vcf_pos > max_bmap_pos)
    if above_count > 0:
        print(f"WARNING: {above_count} VCF positions are above the max B-map start position ({max_bmap_pos}) consider calculating an extended B-map")

    vcf_unique = np.unique(vcf_chroms)
    bmap_unique = np.unique(bmap_chroms)

    missing_in_bmap = set(vcf_unique) - set(bmap_unique)
    missing_in_vcf = set(bmap_unique) - set(vcf_unique)
    if missing_in_bmap:
        print(f"WARNING: The following chromosomes are in the VCF but not in the B-map: {missing_in_bmap}")
    if missing_in_vcf:
        print(f"WARNING: The following chromosomes are in the B-map but not in the VCF: {missing_in_vcf}")

    writer = None
    if args.out:
        out_f = open(args.out, 'w', newline='')
        writer = csv.writer(out_f)
        writer.writerow(['chromosome', 'position', 'B'])

    for chrom in vcf_unique:
        if chrom in missing_in_bmap:
            continue

        mask_v = (vcf_chroms == chrom)
        vcf_pos_chr = vcf_pos[mask_v]

        mask_b = (bmap_chroms == chrom)
        bmap_pos_chr = bmap_pos[mask_b]
        bmap_vals_chr = b_values[mask_b]

        idx = np.searchsorted(bmap_pos_chr, vcf_pos_chr, side='right') - 1
        idx[idx < 0] = 0
        vcf_b_for_chr = bmap_vals_chr[idx]

        if args.out is not None:
            for p, b in zip(vcf_pos_chr, vcf_b_for_chr):
                writer.writerow([chrom, p, b])
        else:
            print(f"{chrom} positions: {vcf_pos_chr}")
            print(f"{chrom} B-values:  {vcf_b_for_chr}")

    if args.out is not None:
        out_f.close()
        print(f"Wrote CSV to {args.out}")
    else: print(f"Skipping save, to save, add --out and --out_binsize")

    print("VCF utilities done")



# 1. Return histogram of B for VCF positions
###     Also give B for each VCF position in file
# 2. Get position of sites with B above or below X, or top N positions of bottom N positions
###     Also recode VCF with those positions if specified


