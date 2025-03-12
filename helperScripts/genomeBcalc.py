from helperScripts.RunBCalcScripts.process_single_chunk import process_single_chunk
from helperScripts.RunBCalcScripts.bedgffHandler import bedgffHandler
from helperScripts.RunBCalcScripts.calculateLPerChunk import calculateLPerChunk
from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.RunBCalcScripts.recmapHandler import recmapHandler
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import os

def genomeBcalc(args):    
    file_path, chr_start, chr_end, calc_start, calc_end, chunk_size, precise_chunks, out, silent = args.bedgff_path, args.chr_start, args.chr_end, args.calc_start, args.calc_end, args.chunk_size, args.precise_chunks, args.out, args.silent

    print(f"= Calculating relative diversity (B) for all neutral sites across the genome. = = =")
    if not args.silent: 
        print(f"====== P A R A M E T E R S =========================")
        print(f"BED/GFF file for regions under selection: {file_path}")
        print(f"First position in chromosome: {calc_start}")
        print(f"Last position in chromosome: {calc_end}")
        print(f"Size of chunks to calculate B in per iteration: {chunk_size}bp")
        print(f"Number of adjacent chunks to calculate B precisely for: {precise_chunks}")

    blockstart, blockend = bedgffHandler(file_path) # Read BED/GFF, return start and end of conserved elements

    if not silent: print(f"====== S T A R T I N G ===== C A L C ===============")

    b_values = np.ones(calc_end + 1 - calc_start, dtype=np.float64) # Initialize array of B values
    for s, e in zip(blockstart, blockend): # Converts gene sites to NaN
        b_values[s - calc_start : e - calc_start + 1] = np.nan

    lperchunk = calculateLPerChunk(chunk_size, blockstart, blockend, chr_start, chr_end) # Cumulative conserved length in each chunk

    if args.rec_map: # Process recombination map if provided
        print(f"Using recombination (crossover) map from {args.rec_map}")
        rec_rate_per_chunk = recmapHandler(args.rec_map, chr_start, chr_end, chunk_size)
    else:
        rec_rate_per_chunk = None

    if args.gc_map:
        print(f"Using gene conversion map from {args.gc_map}")
        gc_rate_per_chunk = recmapHandler(args.gc_map, chr_start, chr_end, chunk_size)
    else:
        gc_rate_per_chunk = None

    if not silent: print(f"====== R E S U L T S == P E R == C H U N K =========")

    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size
    calc_chunk_start = (calc_start - chr_start) // chunk_size
    calc_chunk_end = (calc_end - chr_start) // chunk_size
    calc_chunks = np.arange(calc_chunk_start,calc_chunk_end + 1) # Relevant chunks to calculate B for based on calc_start and calc_end

    with ThreadPoolExecutor() as executor:
        results = [executor.submit(process_single_chunk, chunk_num, 
                                   chunk_size, blockstart, blockend, chr_start, chr_end, calc_start, 
                                   calc_end, num_chunks, precise_chunks, lperchunk, b_values, rec_rate_per_chunk, gc_rate_per_chunk, silent)
            for chunk_num in calc_chunks]
    
    if not silent: 
        print(f"====== F I N I S H E D ===== C A L C ===============")
        print(f"====== R E S U L T S ====== S U M M A R Y ==========")
        print(f"Cumulative length of regions under selection: {int(sum(lperchunk))}bp ({round((sum(lperchunk)/(calc_end - calc_start))*100,2)}%)")
        print(f"Mean B of neutral sites across genome: {b_values[~np.isnan(b_values)].mean()}")
        if args.rec_map: # Process recombination map if provided
            print(f"Calculated using recombination (crossover) map, with rates averaged within {chunk_size}bp chunks")
        if args.gc_map: # Process recombination map if provided
            print(f"Calculated using gene conversion map, with rates averaged within {chunk_size}bp chunks")    

    positions = np.arange(calc_start, calc_end + 1)
    if args.pop_change:
        b_values = get_Bcur(b_values)
        if not silent: print("Demographic change applied to B-calculation")
    output_data = np.column_stack((positions, b_values))

    if args.out is not None:
        csv_file = args.out  # This might be "b_values.csv" or a custom path
        np.savetxt(         # Write to CSV
            csv_file, 
            output_data,
            delimiter=",", 
            header="Position,B", 
            fmt=("%d", "%.6f"),  # first column = integer, second = float w/ 6 decimals
            comments=""
)
        print(f"Saved B values to: {os.path.abspath(csv_file)}")
    else:
        if not args.silent: print("No output CSV requested; skipping save.")

    return output_data