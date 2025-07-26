import numpy as np
import csv
from Bvalcalc.utils.load_vcf import load_vcf
from Bvalcalc.utils.load_Bmap import load_Bmap

def positionsBstats(args, Bmap_path):

    if not args.quiet: print(f"= B-map utilities = = =")
    # Load VCF and B-map
    vcf_chroms, vcf_pos = load_vcf(args.positions)
    print(f"loading B-map: {Bmap_path}")
    bmap_chroms, bmap_pos, b_values = load_Bmap(file_path=Bmap_path)

    # Header for retrieval
    if not args.quiet:
        print("====== R E T R I E V I N G === B - V A L U E S =====")

    # Identify unique chromosomes
    bmap_unique = np.unique(bmap_chroms)
    vcf_unique  = np.unique(vcf_chroms)

    # Report mismatches
    missing_bmap = set(vcf_unique.astype(str)) - set(bmap_unique.astype(str))
    missing_vcf  = set(bmap_unique.astype(str)) - set(vcf_unique.astype(str))
    if missing_bmap:
        print("WARNING: Chromosomes in VCF/CSV but not B-map: " +
              ", ".join(sorted(missing_bmap)))
    if missing_vcf:
        print("WARNING: Chromosomes in B-map but not VCF/CSV: " +
              ", ".join(sorted(missing_vcf)))

    # Filter out VCF entries not in B-map
    if missing_bmap:
        keep_mask = np.isin(vcf_chroms, bmap_unique)
        vcf_chroms, vcf_pos = vcf_chroms[keep_mask], vcf_pos[keep_mask]
        vcf_unique = np.unique(vcf_chroms)

    # Prepare flat arrays
    n = vcf_pos.size
    flat_b     = np.empty(n, dtype=np.float64)
    flat_chrom = np.empty(n, dtype='<U20')

    # Prepare output file if requested
    writer = None
    if args.out:
        out_f = open(args.out, 'w', newline='')
        if not args.bcftools_format:
            # standard CSV with header
            writer = csv.writer(out_f)
            writer.writerow(['chromosome', 'position', 'B'])
        # else: plain "chromosome:position" lines, no header

    # Map B-values per chromosome
    for chrom in vcf_unique:
        mask_v  = (vcf_chroms == chrom)
        pos_chr = vcf_pos[mask_v]
        mask_b  = (bmap_chroms == chrom)
        starts, vals = bmap_pos[mask_b], b_values[mask_b]

        # Overflow warning
        max_start = int(starts.max())
        above_count_chr = int(np.count_nonzero(pos_chr > max_start))
        if above_count_chr > 0:
            print(f"WARNING: {above_count_chr} VCF/CSV positions in {chrom} "
                  f"are above max B-map start ({max_start})")

        # Assign B-values
        idx = np.searchsorted(starts, pos_chr, side='right') - 1
        idx[idx < 0] = 0
        flat_b[mask_v]     = vals[idx]
        flat_chrom[mask_v] = chrom

        # Write out
        if writer:
            for p, b in zip(pos_chr, vals[idx]):
                if (args.out_minimum is None or b >= args.out_minimum) and \
                   (args.out_maximum is None or b <= args.out_maximum):
                    writer.writerow([chrom, p, b])
        elif args.bcftools_format and args.out:
            for p, b in zip(pos_chr, vals[idx]):
                if (args.out_minimum is None or b >= args.out_minimum) and \
                   (args.out_maximum is None or b <= args.out_maximum):
                    out_f.write(f"{chrom}:{p}\n")

    # Summary stats
    filter_mask = np.ones(n, dtype=bool)
    if args.out_minimum is not None:
        filter_mask &= (flat_b >= args.out_minimum)
    if args.out_maximum is not None:
        filter_mask &= (flat_b <= args.out_maximum)

    filtered_b     = flat_b[filter_mask]
    filtered_chrom = flat_chrom[filter_mask]
    filtered_pos   = vcf_pos[filter_mask]
    n_filtered     = filtered_b.size

    if not args.quiet:
        print("====== R E S U L T S ====== S U M M A R Y ==========")
    if n_filtered > 0:
        mean_B    = float(filtered_b.mean())
        idx_max   = int(filtered_b.argmax())
        idx_min   = int(filtered_b.argmin())
        max_B     = float(filtered_b[idx_max])
        min_B     = float(filtered_b[idx_min])
        pos_max   = int(filtered_pos[idx_max])
        pos_min   = int(filtered_pos[idx_min])
        chrom_max = filtered_chrom[idx_max]
        chrom_min = filtered_chrom[idx_min]

        if not args.quiet:
            print(f"Mean B across filtered sites: {mean_B:.6f}")
            print(f"Max B across filtered sites: {max_B:.6f} at {chrom_max}:{pos_max}")
            print(f"Min B across filtered sites: {min_B:.6f} at {chrom_min}:{pos_min}")
    else:
        if not args.quiet:
            print("No B-values to summarize after applying thresholds.")

    # Close file if open
    if args.out:
        out_f.close()
        if not args.quiet:
            print(
                f"Wrote CSV to {args.out}"
                if not args.bcftools_format
                else f"Wrote regions to {args.out}"
            )
    else:
        if not args.quiet:
            print("Skipping save (use --out to write CSV or regions file)")

    return filtered_b, filtered_chrom
