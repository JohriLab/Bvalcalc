from Bvalcalc.core.helpers.process_single_chunk import process_single_chunk
from Bvalcalc.core.helpers.calc_L_per_chunk import calculate_L_per_chunk
from Bvalcalc.core.helpers.demography_helpers import get_Bcur
from Bvalcalc.utils.load_rec_map import load_rec_map
from Bvalcalc.utils.bin_outputs import bin_outputs
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
import numpy as np
import os
import sys

def chromBcalc(args, blockstart, blockend, chromosome, unlinked_B, prior_pos = None, prior_b = None, calc_start=None, calc_end=None, chr_size=None, caller="regionBcalc", write_header=True):    
    #Shared arguments between genomeBcalc and regionBcalc
    file_path, chunk_size, precise_chunks, hri, quiet, verbose = args.bedgff, args.chunk_size, args.precise_chunks, args.hri, args.quiet, args.verbose
    
    num_blocks = len(blockstart)
    # Auto-adjust chunk size for large datasets (only if user hasn't manually set chunk_size)
    if args.chunk_size is None:  # If they did not explicitly provide --chunk_size
        # Set default chunk size
        chunk_size = 20000
        original_chunk_size = chunk_size
        if num_blocks > 250000:
            chunk_size = 1000  # Use 1kb chunks for extremely massive datasets
            if not quiet:
                print(f"Extremely massive dataset detected ({num_blocks} blocks). Auto-adjusting chunk size from {original_chunk_size} to {chunk_size} bp for memory efficiency. Use --chunk_size to override.")
        elif num_blocks > 125000:
            chunk_size = 2000  # Use 2kb chunks for massive datasets
            if not quiet:
                print(f"Massive dataset detected ({num_blocks} blocks). Auto-adjusting chunk size from {original_chunk_size} to {chunk_size} bp for memory efficiency. Use --chunk_size to override.")
        elif num_blocks > 50000:
            chunk_size = 5000  # Use 5kb chunks for very large datasets
            if not quiet:
                print(f"Very large dataset detected ({num_blocks} blocks). Auto-adjusting chunk size from {original_chunk_size} to {chunk_size} bp for memory efficiency. Use --chunk_size to override.")
        elif num_blocks > 25000:
            chunk_size = 10000  # Use 10kb chunks for large datasets
            if not quiet:
                print(f"Large dataset detected ({num_blocks} blocks). Auto-adjusting chunk size from {original_chunk_size} to {chunk_size} bp for memory efficiency. Use --chunk_size to override.")
    elif not quiet and num_blocks > 25000:
        print(f"Large dataset detected ({num_blocks} blocks) but using user-specified chunk size of {chunk_size} bp.")
    #Arguments specific to regionBcalc
    if caller == "regionBcalc":
        calc_start, calc_end = calc_start, calc_end
        if calc_end > blockend[-1]:
            chr_size = calc_end
        else:
            chr_size = None

    if not args.quiet: 
        print(f"====== P A R A M E T E R S =========================")
        print(f"BED/GFF file for regions under selection: {file_path}")
        if chr_size is not None: print(f"Last position in chromosome {chromosome}: {chr_size}")
        print(f"Size of chunks to calculate B in per iteration: {chunk_size}bp")
        print(f"Number of adjacent chunks to calculate B precisely for: {precise_chunks}")

    if chr_size is not None and chr_size < blockend[-1]:
        raise ValueError(f"chr_size provided is less than gene position for chromosome {chromosome}")
    if chr_size is None: # Default chr_size to last value in blockend if not given
        if len(blockend) == 0 and caller != "regionBcalc":
            raise ValueError("chr_size was not provided for chromosome: {chromosome} and gene position ends not computed. Check BED/GFF input, and specify chr_size if needed")
        chr_size = blockend[-1]
        if calc_end is None and not args.quiet:
            print(f"No --chr_size provided for chromosome: {chromosome}. Using last position in BED/GFF: {chr_size}")

    if not quiet: print(f"====== S T A R T I N G ===== C A L C ===============")
    if calc_start is None and calc_end is None:
        if not quiet: print(f"Calculating B for entire chromosome, to only calculate for a subregion, use --calc_start and --calc_end")
    if calc_start is None:
        calc_start = 1
    if calc_end is None:
        calc_end = chr_size

    chr_start = 1 # Currently hardcoded, can change if needed
    num_chunks = (chr_size - chr_start + chunk_size - 1) // chunk_size

    calc_chunk_start = (calc_start - chr_start) // chunk_size
    calc_chunk_end = (calc_end - chr_start) // chunk_size
    calc_chunks = np.arange(calc_chunk_start,calc_chunk_end + 1) # Relevant chunks to calculate B for based on calc_start and calc_end

    b_values = np.ones(chr_size + 2 - chr_start, dtype=np.float64) # Initialize array of B values

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

    import gc
    BATCH_SIZE = args.chunk_batch_size
    total_chunks = len(calc_chunks)
    completed = 0

    for batch_start in range(0, total_chunks, BATCH_SIZE):
        batch = calc_chunks[batch_start : batch_start + BATCH_SIZE]
        with ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(process_single_chunk, chunk_idx,
                                chunk_size, blockstart, blockend, chr_start, chr_size, calc_start,
                                calc_end, num_chunks, precise_chunks, lperchunk, b_values,
                                rec_rate_per_chunk, gc_rate_per_chunk, hri, quiet, verbose, unlinked_B): chunk_idx
                for chunk_idx in batch
            }
            if not quiet and not verbose:
                for future in as_completed(futures):
                    completed += 1
                    progress = int((completed / total_chunks) * 100)
                    sys.stdout.write(f"\rProgress ({chromosome}): {progress}% ({completed}/{total_chunks} chunks [{chunk_size}])")
                    sys.stdout.flush()
                # After batch is done, cleanup
                print()  # Move to the next line after progress printing
                del futures
                gc.collect()

    b_values = b_values[calc_start:(calc_end+1)] # Trim b_values array to only calculated region
    b_values = b_values * unlinked_B
    # Apply prior B values as multipliers after chunk processing
    if prior_pos is not None and prior_b is not None: # Apply prior B values to ranges
        sorted_indices = np.argsort(prior_pos) # Sort positions and B values
        sorted_pos, sorted_b = prior_pos[sorted_indices], prior_b[sorted_indices]
        for i in range(len(sorted_pos)): # Apply B values from each start position to next
            start_pos, b_val = sorted_pos[i], sorted_b[i]
            end_pos = sorted_pos[i + 1] if i < len(sorted_pos) - 1 else calc_end + 1 # Next position or end
            range_start, range_end = max(start_pos, calc_start), min(end_pos, calc_end + 1) # Clip to calc region
            if range_start < range_end:
                idx_start, idx_end = range_start - calc_start, range_end - calc_start # Convert to array indices
                b_values[idx_start:idx_end] *= b_val # Multiply by prior B values

    if hri and rec_rate_per_chunk is not None: # If --hri is active
        from Bvalcalc.core.helpers.extend_hri_regions_correction import extend_hri_regions_correction
        hri_extended_starts, hri_extended_ends = extend_hri_regions_correction(b_values, rec_rate_per_chunk, chunk_size, chr_start, calc_start, calc_end, hri_r_threshold = 0.1) # Extend HRI regions until B > B' to avoid sharp decrease in B at the border between normal and HRI regions. See manuscript.
    else:
        hri_extended_starts, hri_extended_ends = np.array([], dtype=int), np.array([], dtype=int)
    
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
        if prior_pos is not None and prior_b is not None: # Report mean prior B if provided
            mean_prior_b = np.mean(sorted_b)
            print(f"Mean prior B value: {mean_prior_b:.6f}")
        print(f"Cumulative length of chromosome under selection: {int(sum(lperchunk))}bp ({round((sum(lperchunk)/(chr_size - chr_start + 1))*100,2)}%)")
        print(f"B from unlinked sites calculated for chromosome {chromosome}: {unlinked_B}")
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
        names='Chromosome,Start,B',formats='U20,i8,f8')

    if args.out is not None: # Write to CSVs
        print(f"Writing B output to file...")
        from Bvalcalc.utils.write_chrom_B_to_file import write_chrom_B_to_file
        write_chrom_B_to_file(args.out, output_data, quiet, hri_extended_starts, hri_extended_ends, args.out_binsize, calc_end, write_header=write_header, no_header=getattr(args, 'no_header', False))
        print(f"Appended B values to: {os.path.abspath(args.out)}")
    else:
        if not args.quiet:
            print("No output CSV requested; skipping save.")

    if caller == "regionBcalc":
        if rec_rate_per_chunk is not None:
            rec_rate_per_chunk_in_region = rec_rate_per_chunk[calc_start // chunk_size:] # Slice rec_rate_per_chunk from region start onward
        else: rec_rate_per_chunk_in_region = None
        return output_data, block_ranges, rec_rate_per_chunk_in_region, chunk_size
    else: #caller is genomeBcalc
        return