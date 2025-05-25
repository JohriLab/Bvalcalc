from bvalcalc.core.helpers.process_single_chunk import process_single_chunk
from bvalcalc.utils.load_bed_gff import load_bed_gff
from bvalcalc.core.helpers.calc_L_per_chunk import calculate_L_per_chunk
from bvalcalc.core.helpers.demography_helpers import get_Bcur
from bvalcalc.utils.load_rec_map import load_rec_map
from bvalcalc.utils.bin_outputs import bin_outputs
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
import numpy as np
import os
import sys

def chromBcalc(args, blockstart, blockend, chromosome, unlinked_B, prior_pos = None, prior_b = None, calc_start=None, calc_end=None, chr_size=None, caller="regionBcalc"):    
    #Shared arguments between genomeBcalc and regionBcalc
    file_path, chunk_size, precise_chunks, quiet, verbose = args.bedgff_path, args.chunk_size, args.precise_chunks, args.quiet, args.verbose
    #Arguments specific to regionBcalc
    if caller == "regionBcalc":
        calc_start, calc_end = calc_start, calc_end
        chr_size = None

    if not args.quiet: 
        print(f"====== P A R A M E T E R S =========================")
        print(f"BED/GFF file for regions under selection: {file_path}")
        if chr_size is not None: print(f"Last position in chromosome: {calc_end}")
        print(f"Size of chunks to calculate B in per iteration: {chunk_size}bp")
        print(f"Number of adjacent chunks to calculate B precisely for: {precise_chunks}")

    if chr_size is not None and chr_size < blockend[-1]:
        raise ValueError(f"chr_size provided is less than gene position for chromosome {chromosome}")
    if chr_size is None: # Default chr_size to last value in blockend if not given
        if len(blockend) == 0 and caller != "regionBcalc":
            raise ValueError("chr_size was not provided and gene position ends not computed. Check BED/GFF input, and specify chr_size if needed")
        chr_size = blockend[-1]
        if calc_end is None and not args.quiet:
            print(f"No --chr_size provided. Using last position in BED/GFF: {chr_size}")
    if calc_start is None:
        calc_start = 1
    if calc_end is None:
        calc_end = chr_size

    if not quiet: print(f"====== S T A R T I N G ===== C A L C ===============")
    if calc_start is None and calc_end is None:
        if not quiet: print(f"Calculating B for entire chromosome, to only calculate for a subregion, use --calc_start and --calc_end")

    chr_start = 1 # Currently hardcoded, can change if needed
    num_chunks = (chr_size - chr_start + chunk_size - 1) // chunk_size

    calc_chunk_start = (calc_start - chr_start) // chunk_size
    calc_chunk_end = (calc_end - chr_start) // chunk_size
    calc_chunks = np.arange(calc_chunk_start,calc_chunk_end + 1) # Relevant chunks to calculate B for based on calc_start and calc_end

    b_values = np.ones(chr_size + 2 - chr_start, dtype=np.float64) # Initialize array of B values
    if prior_pos is not None and prior_b is not None: # If we have prior B map, overwrite those positions' B values
        idx = np.asarray(prior_pos, dtype=int)
        calc_mask = (idx >= calc_start) & (idx <= calc_end)
        idx = idx[calc_mask] # filter to only those within [calc_start, calc_end]
        bprior = np.asarray(prior_b, dtype=b_values.dtype)[calc_mask]
        b_values[idx] = bprior

    lperchunk = calculate_L_per_chunk(chunk_size, blockstart, blockend, chr_start, chr_size) # Cumulative conserved length in each chunk

    if args.rec_map: # Process recombination map if provided
        if not quiet: print(f"Using recombination (crossover) map from {args.rec_map}")
        rec_rate_per_chunk = load_rec_map(args.rec_map, chr_start, chr_size, chunk_size, chromosome)
    else:
        rec_rate_per_chunk = None

    if args.gc_map:
        if not quiet: print(f"Using gene conversion map from {args.gc_map}")
        gc_rate_per_chunk = load_rec_map(args.gc_map, chr_start, chr_size, chunk_size, chromosome)
    else:
        gc_rate_per_chunk = None

    if verbose: print(f"====== R E S U L T S == P E R == C H U N K =========")
    elif not quiet: print(f"To print per-chunk summaries, add --verbose.")

    with ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(process_single_chunk, chunk_idx,
                            chunk_size, blockstart, blockend, chr_start, chr_size, calc_start,
                            calc_end, num_chunks, precise_chunks, lperchunk, b_values,
                            rec_rate_per_chunk, gc_rate_per_chunk, quiet, verbose): chunk_idx
            for chunk_idx in calc_chunks
        }
        if not quiet and not verbose:
            completed = 0 # Print progress
            for future in as_completed(futures):
                completed += 1
                progress = int((completed / len(calc_chunks)) * 100)
                sys.stdout.write(f"\rProgress: {progress}% ({completed}/{len(calc_chunks)} chunks [{chunk_size}])")
                sys.stdout.flush()
            print()  # Move to the next line after progress printing

    b_values = b_values[calc_start:(calc_end+1)] # Trim b_values array to only calculated region
    b_values = b_values * unlinked_B
    
    if not quiet: 
        print(f"====== F I N I S H E D ===== C A L C ===============")
        print(f"====== R E S U L T S ====== S U M M A R Y ==========")
                # Total genic bases within calc_start to calc_end
        calc_selected_length = 0
        for start, end in zip(blockstart, blockend):
            # Find overlap between this block and the calculated region
            overlap_start = max(start, calc_start)
            overlap_end = min(end, calc_end)
            if overlap_start <= overlap_end:
                calc_selected_length += (overlap_end - overlap_start + 1)
        print(f"Cumulative length of calculated region under selection: {calc_selected_length}bp "f"({round((calc_selected_length / (calc_end - calc_start + 1)) * 100, 2)}%)")
        print(f"Cumulative length of chromosome under selection: {int(sum(lperchunk))}bp ({round((sum(lperchunk)/(chr_size - chr_start + 1))*100,2)}%)")
        print(f"B from unlinked sites for chromosome {chromosome}: {unlinked_B}")
        if caller == "genomeBcalc": print(f"Mean B of neutral sites across chromosome {chromosome}: {b_values[~np.isnan(b_values)].mean()}")
        elif caller == "regionBcalc": print(f"Mean B of neutral sites across specified region: {b_values[~np.isnan(b_values)].mean()}")
        if args.rec_map: # Process recombination map if provided
            print(f"Calculated using recombination (crossover) map, with rates averaged within {chunk_size}bp chunks")
        if args.gc_map: # Process recombination map if provided
            print(f"Calculated using gene conversion map, with rates averaged within {chunk_size}bp chunks")    

    block_ranges = np.column_stack((np.repeat(chromosome, blockstart.shape[0]), blockstart, blockend))

    positions = np.arange(calc_start, calc_end + 1)
    conserved = np.full_like(positions, "N", dtype="<U1")
    for start, end in zip(blockstart, blockend): # Mark conserved regions
        conserved[max(start, calc_start) - calc_start : min(end, calc_end) - calc_start + 1] = "C"

    if args.pop_change:
        b_values = get_Bcur(b_values)
        if not quiet: print("Demographic change applied to B-calculation")

    binned_b_values, binned_positions = bin_outputs(b_values, positions, args.out_binsize)
    chrom_col = np.full(binned_positions.shape, chromosome, dtype="<U20")

    output_data = np.core.records.fromarrays(
    [chrom_col,binned_positions.astype(int),binned_b_values.astype(float)],
    names='Chromosome,Position,B',formats='U20,i8,f8')

    if args.out is not None: # Write to CSV
        print(f"Writing B output to file...")
        with open(args.out, 'a') as f:
            np.savetxt(f, output_data, delimiter=",", fmt="%s,%d,%.6f", comments="")
        print(f"Appended B values to: {os.path.abspath(args.out)}")
    else:
        if not args.quiet:
            print("No output CSV requested; skipping save.")

    if caller == "regionBcalc":
        return output_data, block_ranges
    else: #caller should be genomeBcalc
        return