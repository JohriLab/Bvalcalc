#!/usr/bin/env python3

# ./Bcalc_multiElements.py --chr_start 1 --chr_end 25254535 --file_path ../exampleData/dmel6_2R_genes.csv

#calculate_B_analytically_Eq3_mine_demography_multiple_elements
#This script is to get B values across a genomic element with multiple functional elements.
#Currently the output is in terms of an average of a sliding window
import sys
import math
import numpy as np
from numpy.lib import recfunctions
import csv
from Bcalc_function import calculate_B
from calcLfromB_function import find_minimum_distance_binary
from constants import g, tract_len, r, u, Ncur, Nanc, gamma_cutoff, h, t0, t1, t1half, t2, t3, t4, f0, f1, f2, f3
import argparse

#CLI handling
def main():
    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    parser.add_argument('--chr_start', type=int, required=True, help="Start of chromosome range.")
    parser.add_argument('--chr_end', type=int, required=True, help="End of chromosome range.")
    parser.add_argument('--chunk_size', type=int, default=100000, help="Size of chunks calculated simulataneously (bp), decrease if running out of memory, increase for performance. [100000]")
    parser.add_argument('--flank_len', type=int, default=100000, help="Length of region adjacent to conserved element for which B is calculated (bp). [20000]")
    parser.add_argument('--file_path', type=str, required=True, help="Path to input BED or GFF3 file with conserved regions (e.g. genes).")

    args = parser.parse_args()
 
    runBcalc(args)

#Main function
def runBcalc(args):
    #Define variables and constants:
    file_path = args.file_path #../exampleData/dmel6_2R_genes.csv #../exampleData/test.csv
    chr_start = args.chr_start # 1
    chr_end = args.chr_end # 25254535
    flank_len = args.flank_len # 20000
    chunk_size = args.chunk_size # 100000
    out_folder="example_out"

    ## 1. PARSE INPUT OF BED CONSERVED REGIONS

    blockstart = []
    blockend = []
    seen_blocks = set()

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 2:
                start, end = int(row[1]), int(row[2])
                if (start, end) not in seen_blocks:  # Check if the pair is unique
                    seen_blocks.add((start, end))  # Mark this pair as seen
                    blockstart.append(start)
                    blockend.append(end)
    blockstart = np.array(blockstart)
    blockend = np.array(blockend)
    lengths = blockend - blockstart #XX Need to extend for reverse orientation


# Calculate distances for each length
    flank_distances = np.zeros_like(lengths, dtype=np.int32)
    flank_blockstart = np.zeros_like(blockstart, dtype=np.int32)
    flank_blockend = np.zeros_like(blockend, dtype=np.int32)
    for i, length in enumerate(lengths):
        flank_distances[i] = find_minimum_distance_binary(0.998, length)
        flank_blockstart[i] = blockstart[i] - flank_distances[i]
        flank_blockend[i] = blockend[i] + flank_distances[i]

    print(blockend, flank_blockend)

    # print("Distances for B > 0.998:", distances)


    # 2. CREATE ARRAY OF ALL NEUTRAL SITES (pos, B)
    # Initialize the array for B values (all initially set to 1.0)
    b_values = np.ones(chr_end - chr_start, dtype=np.float64)

    # 3a. EXTRACT POSITIONS 5KB DOWNSTREAM OF ELEMENT
    # Divide positions into chunks

    def generate_chunks(array, chunk_size):
        """Yields chunks of the array one at a time."""
        for i in range(0, len(array), chunk_size):
            yield array[i:i + chunk_size]

    # Generate chunks lazily
    pos_chunk_generator = generate_chunks(np.arange(chr_start, chr_end), chunk_size)

    # Iterate over chunks, calculating B for all neutral sites
    for pos_chunk in pos_chunk_generator:

        # Filter blockstart and blockend to only include where flanking region overlaps with chunk
        relevant_blockregion = (pos_chunk.max() >= flank_blockstart) & (pos_chunk.min() <= flank_blockend)
        relevant_blockstart = blockstart[relevant_blockregion]
        relevant_blockend = blockend[relevant_blockregion]

        # Filter lengths to match the relevant blocks
        relevant_lengths = lengths[relevant_blockregion]

        # Calculate distances and masks for this chunk
        distances_upstream = relevant_blockstart[:, None] - pos_chunk[None, :] 
        distances_downstream = pos_chunk[None, :] - relevant_blockend[:, None] 

        # Masks for flanking sites
        upstream_mask = (pos_chunk < relevant_blockstart[:, None]) & \
                        (pos_chunk > (relevant_blockstart[:, None] - flank_len))
        downstream_mask = (pos_chunk > relevant_blockend[:, None]) & \
                        (pos_chunk < (relevant_blockend[:, None] + flank_len))
        flanking_mask = upstream_mask | downstream_mask
        
        # Combine the distances into a single array
        distances = np.where(flanking_mask, 
                            np.where(upstream_mask, distances_upstream, distances_downstream), 
                            np.nan)
        # Flatten distances and flanking_mask to match the selected elements
        flat_distances = distances[flanking_mask]  # Select distances where mask is True
        flat_lengths = np.repeat(relevant_lengths, flanking_mask.sum(axis=1))[:len(flat_distances)]  # Repeat each length to match the mask, handling edge cases where lengths overshoot
        # Calculate B for the flattened data
        flank_B = calculate_B(flat_distances, flat_lengths)
        # Flatten flanking_mask to get indices of True values
        true_indices = np.where(flanking_mask)
        # Find unique indices and map each to its position in the unique list
        unique_indices, inverse_indices = np.unique(true_indices[1], return_inverse=True)
        # Aggregate B values for unique indices
        aggregated_B = np.ones_like(unique_indices, dtype=np.float64)
        np.multiply.at(aggregated_B, inverse_indices, flank_B)
        # Find global indices for the current chunk
        global_indices = pos_chunk[unique_indices] - chr_start  # Convert chunk indices to global indices
        # Incrementally update the `B` values in the b_values array
        # Memory-efficient multiplication to update `B` in growing b_values array
        for idx, agg_B in zip(global_indices, aggregated_B):
            b_values[idx] *= agg_B

        print(f"Processing chunk: {pos_chunk.min()} - {pos_chunk.max()}")
        print(f"Number of relevant genes: {len(relevant_blockstart)}")
        print(f"Relevant blocks: {relevant_blockstart}, {relevant_blockend}")
        print(f"Aggregated B values for chunk: {aggregated_B}")


if __name__ == "__main__":
    main()

sys.exit()

# Converts B at gene sites to 'nan'
for start, end in zip(blockstart, blockend):
    sites['B'][start - sites['pos'][0]:end - sites['pos'][0]] = np.nan 

nogene_sites = np.sort(sites[~np.isnan(sites['B'])], order=['B']) # Removes gene sites
print(nogene_sites)
sys.exit()

#calculate an average B value over the window with coordinates win_start - win_end
def calculate_Banc_window(win_start, win_end):
    i=win_start
    B_sum = 0.0
    s_tot = 0
    while i <= win_end:
        #!!!if position i is a neutral site: #!!! Check if this site is neutral. Calculate B if it is neutral. Otherwise do not compute B.
            #calculating B at site i from all genomic elements:
        #!!!    B_i = 1.0
        #!!!    For each functional element in the genome: #!!!
        #!!!        B_i = B_i * float(calculate_B(distance_to_element[i], length_of_element[i])) #!!! give the nearest distance from position i to the given element and the length of that element
       #!!!     B_sum = B_sum + B_i
        #!!!    s_tot = s_tot + 1
        i = i + 1
    if s_tot > 0: 
        return(B_sum/s_tot)
    else: #If all sites in this window are selected, then return NA.
        return ("NA")

#Gets the value of B in the current population as a function of B calcualted in the ancestral population (assumed to be in equilibrium).
#As in, this is where we account for a simple single-size change in N
def get_Bcur(Banc):
    R = float(Nanc)/float(Ncur)
    Bcur = (Banc*(1.0 + (R-1.0)*math.exp((-1.0*time_of_change)/Banc))) / (1 + (R-1.0)*math.exp(-1.0*time_of_change))
    return (Bcur)

                                                                                
#Getting an average pi over a window:
result = open("/work/users/p/j/pjohri/BGSCalculator/mytheory/" + out_folder + "/Bvalues_" + DFE + "_" + TIME + "_windowsize_" + str(s_window_size) + ".txt", 'w+')
result.write("start" + '\t' + "end" + '\t' + "B" + '\n')
posn_start = 1
while posn_start <= chr_end:
    result.write(str(posn_start) + '\t' + str(posn_start+s_window_size) + '\t' + str(get_Bcur(calculate_Banc_window(posn_start, posn_start+s_window_size))) + '\n')
    posn_start = posn_start + s_window_size
result.close()

print("done")




# # Extracting neutral sites
# neu_positions = np.array([
#     site for site in sites['pos']
#     if not any(start <= site <= end for start, end in zip(blockstart, blockend))
# ])

# neu_sites = np.zeros(len(neu_positions), dtype=[('pos', 'int32'), ('B', 'f4')])
# neu_sites["pos"] = neu_positions
# neu_sites["B"] = 1