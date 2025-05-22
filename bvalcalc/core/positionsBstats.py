from bvalcalc.utils.load_vcf import load_vcf
from bvalcalc.utils.load_Bmap import load_Bmap
import numpy as np
import csv
import sys

def positionsBstats(args, Bmap_path):    
    # Load info using utils
    vcf_chroms, vcf_pos = load_vcf(args.positions)
    bmap_chroms, bmap_pos, b_values = load_Bmap(file_path=Bmap_path)

    # Check for positions beyond B-map range
    max_bmap_pos = bmap_pos[-1]
    above_count = np.sum(vcf_pos > max_bmap_pos)
    if above_count > 0:
        print(f"WARNING: {above_count} VCF positions are above the max B-map start position ({max_bmap_pos}) consider calculating an extended B-map")

    # Identify unique chromosomes and mismatches
    vcf_unique = np.unique(vcf_chroms)
    bmap_unique = np.unique(bmap_chroms)

    missing_in_bmap = set(vcf_unique) - set(bmap_unique)
    missing_in_vcf = set(bmap_unique) - set(vcf_unique)
    if missing_in_bmap:
        print(f"WARNING: The following chromosomes are in the VCF but not in the B-map: {missing_in_bmap}")
    if missing_in_vcf:
        print(f"WARNING: The following chromosomes are in the B-map but not in the VCF: {missing_in_vcf}")
    print(f"====== R E T R I E V I N G === B - V A L U E S =====")

    # Prepare CSV writer if output requested
    writer = None
    if args.out:
        out_f = open(args.out, 'w', newline='')
        writer = csv.writer(out_f)
        writer.writerow(['chromosome', 'position', 'B'])

    # Collect all B-values, positions, and chroms for summary stats
    all_b_list = []
    all_pos_list = []
    all_chrom_list = []

    # Loop through each chromosome in VCF
    for chrom in vcf_unique:
        if chrom in missing_in_bmap:
            continue

        # Extract VCF positions for this chrom
        mask_v = (vcf_chroms == chrom)
        vcf_pos_chr = vcf_pos[mask_v]

        # Extract B-map starts & values for this chrom
        mask_b = (bmap_chroms == chrom)
        bmap_pos_chr = bmap_pos[mask_b]
        bmap_vals_chr = b_values[mask_b]

        # Map positions to B-values via searchsorted
        idx = np.searchsorted(bmap_pos_chr, vcf_pos_chr, side='right') - 1
        idx[idx < 0] = 0
        vcf_b_for_chr = bmap_vals_chr[idx]

        # Append to summary lists
        all_b_list.append(vcf_b_for_chr)
        all_pos_list.append(vcf_pos_chr)
        all_chrom_list.append(np.array([chrom] * len(vcf_pos_chr), dtype='<U20'))

        # Write or print per-chrom data
        if writer:
            for p, b in zip(vcf_pos_chr, vcf_b_for_chr):
                writer.writerow([chrom, p, b])
        else:
            print(f"{chrom} positions: {vcf_pos_chr}")
            print(f"{chrom} B-values:  {vcf_b_for_chr}")

    # Summary statistics
    print(f"====== R E S U L T S ====== S U M M A R Y ==========")
    if all_b_list:
        all_b = np.concatenate(all_b_list)
        all_pos = np.concatenate(all_pos_list)
        all_chrom = np.concatenate(all_chrom_list)
        mean_B = np.mean(all_b)
        max_B = np.max(all_b)
        min_B = np.min(all_b)
        # find corresponding positions and chroms
        idx_max = np.argmax(all_b)
        idx_min = np.argmin(all_b)
        pos_max = all_pos[idx_max]
        pos_min = all_pos[idx_min]
        chrom_max = all_chrom[idx_max]
        chrom_min = all_chrom[idx_min]
        print(f"Mean B across VCF: {mean_B:.6f}")
        print(f"Max B across VCF: {max_B:.6f} at {chrom_max}:{pos_max}")
        print(f"Min B across VCF: {min_B:.6f} at {chrom_min}:{pos_min}")
    else:
        print("No B-values to summarize.")

    # Close CSV file if used
    if writer:
        out_f.close()
        print(f"Wrote CSV to {args.out}")
    else:
        print(f"Skipping save, to save, add --out and --out_binsize")

    return
