import sys
import math
import numpy as np
from numpy.lib import recfunctions
import csv
from helperScripts.calculate_B import calculate_B
from BvalueCalculator.pythonScripts.helperScripts.findFlankLen import find_minimum_distance_binary
from constants import g, tract_len, r, u, Ncur, Nanc, gamma_cutoff, h, t0, t1, t1half, t2, t3, t4, f0, f1, f2, f3
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import time

#CLI handling
def main():

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    parser.add_argument('--chr_start', type=int, required=True, help="Start of chromosome range.")
    parser.add_argument('--chr_end', type=int, required=True, help="End of chromosome range.")
    parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of chunks calculated simulataneously (bp), decrease if running out of memory, increase for performance. [100000]")
    parser.add_argument('--flank_len', type=int, default=100000, help="Length of region adjacent to conserved element for which B is calculated (bp). [20000]")
    parser.add_argument('--file_path', type=str, required=True, help="Path to input BED or GFF3 file with conserved regions (e.g. genes).")

    args = parser.parse_args()
 
    runBcalc(args)

            # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script completed in {elapsed_time:.2f} seconds.")

    # for pos_chunk in pos_chunk_generator:
# Modify process_single_chunk to accept a single argument (a tuple of all required inputs)
def process_single_chunk(args):
    # Unpack the tuple of arguments
    pos_chunk, flank_blockstart, flank_blockend, blockstart, blockend, lengths, chr_start, b_values = args

    # The rest of the function remains unchanged
    relevant_blockregion = (pos_chunk.max() >= flank_blockstart) & (pos_chunk.min() <= flank_blockend)
    relevant_blockstart = blockstart[relevant_blockregion]
    relevant_blockend = blockend[relevant_blockregion]
    relevant_flank_blockstart = flank_blockstart[relevant_blockregion]
    relevant_flank_blockend = flank_blockend[relevant_blockregion]

    relevant_lengths = lengths[relevant_blockregion]
    distances_upstream = relevant_blockstart[:, None] - pos_chunk[None, :]
    distances_downstream = pos_chunk[None, :] - relevant_blockend[:, None]

    upstream_mask = (pos_chunk < relevant_blockstart[:, None]) & (pos_chunk > (relevant_flank_blockstart[:, None]))
    downstream_mask = (pos_chunk > relevant_blockend[:, None]) & (pos_chunk < (relevant_flank_blockend[:, None]))
    flanking_mask = upstream_mask | downstream_mask

    distances = np.where(
        flanking_mask,
        np.where(upstream_mask, distances_upstream, distances_downstream),
        np.nan,
    )
    flat_distances = distances[flanking_mask]
    flat_lengths = np.repeat(relevant_lengths, flanking_mask.sum(axis=1))[: len(flat_distances)]
    flank_B = calculate_B(flat_distances, flat_lengths)

    true_indices = np.where(flanking_mask)
    unique_indices, inverse_indices = np.unique(true_indices[1], return_inverse=True)
    aggregated_B = np.ones_like(unique_indices, dtype=np.float64)
    np.multiply.at(aggregated_B, inverse_indices, flank_B)

    global_indices = pos_chunk[unique_indices] - chr_start
    for idx, agg_B in zip(global_indices, aggregated_B):
        b_values[idx] *= agg_B

    print(f"Processing chunk: {pos_chunk.min()} - {pos_chunk.max()}")
    return b_values

def generate_chunks(array, chunk_size):
    """Yields chunks of the array one at a time."""
    for i in range(0, len(array), chunk_size):
        yield array[i:i + chunk_size]

#Main function
# Parallel processing of the first four chunks
def runBcalc(args):
    # Same setup as before
    file_path = args.file_path
    chr_start = args.chr_start
    chr_end = args.chr_end
    chunk_size = args.chunk_size

    blockstart, blockend = [], []
    seen_blocks = set()
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 2:
                start, end = int(row[1]), int(row[2])
                if (start, end) not in seen_blocks:
                    seen_blocks.add((start, end))
                    blockstart.append(start)
                    blockend.append(end)

    blockstart = np.array(blockstart)
    blockend = np.array(blockend)
    lengths = blockend - blockstart

    flank_distances = np.zeros_like(lengths, dtype=np.int32)
    flank_blockstart = np.zeros_like(blockstart, dtype=np.int32)
    flank_blockend = np.zeros_like(blockend, dtype=np.int32)
    for i, length in enumerate(lengths):
        flank_distances[i] = find_minimum_distance_binary(0.998, length)
        flank_blockstart[i] = blockstart[i] - flank_distances[i]
        flank_blockend[i] = blockend[i] + flank_distances[i]

    b_values = np.ones(chr_end - chr_start, dtype=np.float64)

    pos_chunks = list(generate_chunks(np.arange(chr_start, chr_end), chunk_size))

    # Prepare arguments as tuples for the function
    chunk_args = [
        (chunk, flank_blockstart, flank_blockend, blockstart, blockend, lengths, chr_start, b_values.copy())
        for chunk in pos_chunks
    ]

    # Use a ProcessPoolExecutor with 8 workers
    with ProcessPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(process_single_chunk, args): idx for idx, args in enumerate(chunk_args)}

        for future in as_completed(futures):
            idx = futures[future]
            try:
                result = future.result()
                print(f"Chunk {idx + 1}: Min B value = {np.min(result)}")
            except Exception as e:
                print(f"Chunk {idx + 1} processing failed: {e}")



if __name__ == "__main__":
    main()