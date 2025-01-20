#!/usr/bin/env python3

# ./Bcalc_multiElements.py --chr_start 1 --chr_end 25254535 --file_path ../exampleData/dmel6_2R_genes.csv

#calculate_B_analytically_Eq3_mine_demography_multiple_elements
#This script is to get B values across a genomic element with multiple functional elements.
#Currently the output is in terms of an average of a sliding window
import sys
from helperScripts.parseArgs import parseArgs
from helperScripts.runBcalc import runBcalc
import time

#CLI handling
def main():
    start_time = time.time()

    args = parseArgs() 
    runBcalc(args)

    print(f"Script completed in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()

sys.exit()


# #calculate an average B value over the window with coordinates win_start - win_end
# def calculate_Banc_window(win_start, win_end):
#     i=win_start
#     B_sum = 0.0
#     s_tot = 0
#     while i <= win_end:
#         #!!!if position i is a neutral site: #!!! Check if this site is neutral. Calculate B if it is neutral. Otherwise do not compute B.
#             #calculating B at site i from all genomic elements:
#         #!!!    B_i = 1.0
#         #!!!    For each functional element in the genome: #!!!
#         #!!!        B_i = B_i * float(calculate_B(distance_to_element[i], length_of_element[i])) #!!! give the nearest distance from position i to the given element and the length of that element
#        #!!!     B_sum = B_sum + B_i
#         #!!!    s_tot = s_tot + 1
#         i = i + 1
#     if s_tot > 0: 
#         return(B_sum/s_tot)
#     else: #If all sites in this window are selected, then return NA.
#         return ("NA")

                                                                                
# #Getting an average pi over a window:
# result = open("/work/users/p/j/pjohri/BGSCalculator/mytheory/" + out_folder + "/Bvalues_" + DFE + "_" + TIME + "_windowsize_" + str(s_window_size) + ".txt", 'w+')
# result.write("start" + '\t' + "end" + '\t' + "B" + '\n')
# posn_start = 1
# while posn_start <= chr_end:
#     result.write(str(posn_start) + '\t' + str(posn_start+s_window_size) + '\t' + str(get_Bcur(calculate_Banc_window(posn_start, posn_start+s_window_size))) + '\n')
#     posn_start = posn_start + s_window_size
# result.close()