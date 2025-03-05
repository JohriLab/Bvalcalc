from helperScripts.RunBCalcScripts.process_single_chunk import process_single_chunk
from helperScripts.RunBCalcScripts.bedgffHandler import bedgffHandler
from helperScripts.RunBCalcScripts.calculateLPerChunk import calculateLPerChunk
from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.RunBCalcScripts.recmapHandler import recmapHandler
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import os

def genomeBcalc(args):    
    file_path, actual_chrstart, actual_chrend, calc_start, calc_end, chunk_size, precise_chunks, out = args.bedgff_path, args.calc_start, args.calc_end, args.chunk_size, args.precise_chunks, args.out

    print(f"Calculating relative diversity (B) for all neutral sites across the genome...")
    print(f"====== P A R A M E T E R S =========================")
    print(f"BED/GFF file for regions under selection: {file_path}")
    print(f"First position in chromosome: {calc_start}")
    print(f"Last position in chromosome: {calc_end}")
    print(f"Size of chunks to calculate B in per iteration: {chunk_size}bp")
    print(f"Number of adjacent chunks to calculate B precisely for: {precise_chunks}")

    blockstart, blockend = bedgffHandler(file_path) # Read BED/GFF, return start and end of conserved elements
    b_values = np.ones(calc_end - calc_start, dtype=np.float64) # Initialize array of B values
    for s, e in zip(blockstart, blockend): # Converts gene sites to NaN
        b_values[s - calc_start : e - calc_start + 1] = np.nan

    lperchunk = calculateLPerChunk(chunk_size, blockstart, blockend, calc_start, calc_end) # Cumulative conserved length in each chunk

    if args.rec_map: # Process recombination map if provided
        print(f"Using recombination map from {args.rec_map}")
        rec_rate_per_chunk = recmapHandler(args.rec_map, calc_start, calc_end, chunk_size)
    else:
        rec_rate_per_chunk = None

    print(f"====== R E S U L T S == P E R == C H U N K =========")

    num_chunks = (calc_end - calc_start + chunk_size - 1) // chunk_size
    with ThreadPoolExecutor() as executor:
        results = [executor.submit(process_single_chunk, chunk_num, 
                                   chunk_size, blockstart, blockend, calc_start, 
                                   calc_end, num_chunks, precise_chunks, lperchunk, b_values, rec_rate_per_chunk)
            for chunk_num in range(num_chunks)]
    
    print(f"====== R E S U L T S ====== S U M M A R Y ==========")
    print(f"Mean B of neutral sites across genome: {b_values[~np.isnan(b_values)].mean()}")
    print(f"Cumulative length of regions under selection: {int(sum(lperchunk))}bp ({round((sum(lperchunk)/(calc_end - calc_start))*100,2)}%)")

    positions = np.arange(calc_start, calc_end)
    output_data = np.column_stack((positions, b_values))

    if args.out is not None:
        csv_file = args.out  # This might be "b_values.csv" or a custom path

        # Combine positions and b_values into two columns

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
        return output_data