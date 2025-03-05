from helperScripts.RunBCalcScripts.process_single_chunk import process_single_chunk
from helperScripts.RunBCalcScripts.bedgffHandler import bedgffHandler
from helperScripts.RunBCalcScripts.calculateLPerChunk import calculateLPerChunk
from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.RunBCalcScripts.recmapHandler import recmapHandler
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import os

def genomeBcalc(args):    
    file_path, chr_start, chr_end, chunk_size, precise_chunks, out = args.bedgff_path, args.chr_start, args.chr_end, args.chunk_size, args.precise_chunks, args.out

    print(f"Calculating relative diversity (B) for all neutral sites across the genome...")
    print(f"====== P A R A M E T E R S =========================")
    print(f"BED/GFF file for regions under selection: {file_path}")
    print(f"First position in chromosome: {chr_start}")
    print(f"Last position in chromosome: {chr_end}")
    print(f"Size of chunks to calculate B in per iteration: {chunk_size}bp")
    print(f"Number of adjacent chunks to calculate B precisely for: {precise_chunks}")

    blockstart, blockend = bedgffHandler(file_path) # Read BED/GFF, return start and end of conserved elements
    b_values = np.ones(chr_end - chr_start, dtype=np.float64) # Initialize array of B values
    for s, e in zip(blockstart, blockend): # Converts gene sites to NaN
        b_values[s - chr_start : e - chr_start + 1] = np.nan

    lperchunk = calculateLPerChunk(chunk_size, blockstart, blockend, chr_start, chr_end) # Cumulative conserved length in each chunk

    if args.rec_map: # Process recombination map if provided
        print(f"Using recombination map from {args.rec_map}")
        rec_rate_per_chunk = recmapHandler(args.rec_map, chr_start, chr_end, chunk_size)
    else:
        rec_rate_per_chunk = None

    print(f"====== R E S U L T S == P E R == C H U N K =========")

    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size
    with ThreadPoolExecutor() as executor:
        results = [executor.submit(process_single_chunk, x, chunk_size, blockstart, blockend, chr_start, 
                                   chr_end, num_chunks, precise_chunks, lperchunk, b_values, rec_rate_per_chunk)
            for x in range(num_chunks)]
    
    print(f"====== R E S U L T S ====== S U M M A R Y ==========")
    print(f"Mean B of neutral sites across genome: {b_values[~np.isnan(b_values)].mean()}")
    print(f"Cumulative length of regions under selection: {int(sum(lperchunk))}bp ({round((sum(lperchunk)/(chr_end - chr_start))*100,2)}%)")

    if args.out is not None:
        csv_file = args.out  # This might be "b_values.csv" or a custom path

        # Combine positions and b_values into two columns
        positions = np.arange(1, 200000)
        output_data = np.column_stack((positions, b_values))

        # Write to CSV
        np.savetxt(
            csv_file, 
            output_data,
            delimiter=",", 
            header="Position,B", 
            fmt=("%d", "%.6f"),  # first column = integer, second = float w/ 6 decimals
            comments=""
)
        print(f"Saved B values to: {os.path.abspath(csv_file)}")
    else:
        print("No output CSV requested; skipping save.")

    if args.pop_change:
        return get_Bcur(b_values)
    else:
        return b_values