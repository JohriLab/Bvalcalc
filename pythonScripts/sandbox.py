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
from helperScripts.findFlankLen import findFlankLen
from helperScripts.process_single_chunk import process_single_chunk
from helperScripts.bedgffHandler import bedgffHandler
from constants import g, tract_len, r, u, Ncur, Nanc, gamma_cutoff, h, t0, t1, t1half, t2, t3, t4, f0, f1, f2, f3
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing as mp
import time

#CLI handling
def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    parser.add_argument('--chr_start', type=int, required=True, help="Start of chromosome range.")
    parser.add_argument('--chr_end', type=int, required=True, help="End of chromosome range.")
    parser.add_argument('--chunk_size', type=int, default=100000, help="Size of chunks calculated simulataneously (bp), decrease if running out of memory, increase for performance. [100000]")
    parser.add_argument('--flank_len', type=int, default=100000, help="Length of region adjacent to conserved element for which B is calculated (bp). [20000]")
    parser.add_argument('--file_path', type=str, required=True, help="Path to input BED or GFF3 file with conserved regions (e.g. genes).")

    args = parser.parse_args()
 
    runBcalc(args)

    end_time = time.time()
    print(f"Script completed in {end_time - start_time:.2f} seconds.")

#Main function
def runBcalc(args):    
    file_path, chr_start, chr_end, chunk_size = args.file_path, args.chr_start, args.chr_end, args.chunk_size

    # Read BED/GFF and return genes and relevant flanking regions for calculating B
    blockstart, blockend, lengths, flank_blockstart, flank_blockend =  \
        bedgffHandler(file_path) 

    # Initialize the array for B values (all initially set to 1.0)
    b_values = np.ones(chr_end - chr_start, dtype=np.float64)
    
    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size

    with ThreadPoolExecutor() as executor:
        results = [executor.submit(process_single_chunk, x, chunk_size, flank_blockstart, flank_blockend,
                blockstart, blockend, lengths, chr_start, chr_end, b_values)
                for x in range(num_chunks)]

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