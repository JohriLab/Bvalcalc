#calculate_B_analytically_Eq3_mine_demography_multiple_elements
#This script is to get B values across a genomic element with multiple functional elements.
#Currently the output is in terms of an average of a sliding window

import sys
import math
import numpy as np
from numpy.lib import recfunctions
import csv
import timeit
from Bcalc_function import calculate_B, vectorized_B
from constants import g, tract_len, r, u, l, U, Ncur, Nanc, gamma_cutoff, h, t0, t1, t1half, t2, t3, t4, f0, f1, f2, f3


#Define variables and constants:
out_folder="droso_single_exon_gc_10kb_decline10x"
chr_start = 900
chr_end = 10000
flank_len = 10

#Parameters of genome architecture:
# !!!Read in a bed file with positions of functional elements and the last position (i.e., where the genome ends)#!!!
# !!!Save in memory, the positions of the selected sites. #!!! No need to save the positions of the neutral sites (might take upp too much memory).
# !!!last_position = 10000.0 #(*Full length of the chromosome; the last position. Get this from the user or the input file.!!!*)

#Parameters of instantaneous change in demographic history:

TIME="T_1" #T_0_1/T_0_5/T_1
time_of_change=1.0 #0.1/0.5/1(This is the time of change in 2Ncur generations in the past.)

#Parameters of the DFE:
DFE="DFE3"
s_window_size = 100

#Constants that we do not need from the user:
pi_anc = 4*Nanc*u #(*Expected nucleotide diversity under neutrality*)
blockstart = []
blockend = []

## 1. PARSE INPUT OF BED CONSERVED REGIONS
with open('../exampleData/test.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if len(row) >= 2:
            blockstart.append(int(row[1]))
            blockend.append(int(row[2]))

blockstart = np.array(blockstart)
blockend = np.array(blockend)

# 2. CREATE 2D ARRAY OF ALL NEUTRAL SITES (pos, B)
sites = np.zeros(chr_end - chr_start, dtype=[('pos', 'int32'), ('B', 'f4')])
sites['pos'] = np.arange(chr_start, chr_end)
sites["B"] = 1

# 3a. EXTRACT POSITIONS 5KB DOWNSTREAM OF ELEMENT

# Convert blockstart, blockend, and sites['pos'] to NumPy arrays
pos_array = sites['pos']
blockstart_array = np.array(blockstart)
blockend_array = np.array(blockend)

# Broadcast blockstart and blockend arrays to match the shape of pos_array
lengths = blockend_array - blockstart_array
# Broadcast blockstart and blockend arrays to match the shape of pos_array
distances_upstream = blockstart_array[:, None] - pos_array[None, :]
distances_downstream = pos_array[None, :] - blockend_array[:, None]
# Combine the distances into a single array
# Masks for flanking sites
upstream_mask = (pos_array < blockstart_array[:, None]) & (pos_array > (blockstart_array[:, None] - flank_len))
downstream_mask = (pos_array > blockend_array[:, None]) & (pos_array < (blockend_array[:, None] + flank_len))

# Combine the upstream and downstream masks
flanking_mask = upstream_mask | downstream_mask

# Combine the distances into a single array
distances = np.where(flanking_mask, 
                     np.where(upstream_mask, distances_upstream, distances_downstream), 
                     0)

# Flatten distances and flanking_mask to match the selected elements
flat_distances = distances[flanking_mask]  # Select distances where mask is True
flat_lengths = np.repeat(lengths, flanking_mask.sum(axis=1))  # Repeat each length to match the mask
flat_lengths = flat_lengths[:len(flat_distances)]  # Handle edge cases where lengths might overshoot

# Calculate B for the flattened data
flank_B = calculate_B(flat_distances, flat_lengths)
# Flatten flanking_mask to get indices of True values
true_indices = np.where(flanking_mask)

# Find unique indices and map each to its position in the unique list
unique_indices, inverse_indices = np.unique(true_indices[1], return_inverse=True)

aggregated_B = np.ones_like(unique_indices, dtype=np.float64)

np.multiply.at(aggregated_B, inverse_indices, flank_B)

# Update 'B' values in the original NumPy array
sites['B'][unique_indices] = aggregated_B

# Converts B at gene sites to 'nan'
for start, end in zip(blockstart, blockend):
    sites['B'][start - sites['pos'][0]:end - sites['pos'][0]] = np.nan 

nogene_sites = np.sort(sites[~np.isnan(sites['B'])], order=['B']) # Removes gene sites
# print(nogene_sites)
# print(sites)
print(calculate_B(2,20000))
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
