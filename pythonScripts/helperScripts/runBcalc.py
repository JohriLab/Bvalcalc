from helperScripts.process_single_chunk import process_single_chunk
from helperScripts.bedgffHandler import bedgffHandler
from helperScripts.calculateLPerChunk import calculateLPerChunk
from helperScripts.demographyHelpers import get_Bcur
from constants import g, tract_len, r, u, Ncur, Nanc, gamma_cutoff, h, t0, t1, t1half, t2, t3, t4, f0, f1, f2, f3, time_of_change
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

#Main function
def runBcalc(args):    
    file_path, chr_start, chr_end, chunk_size, precise_chunks = args.file_path, args.chr_start, args.chr_end, args.chunk_size, args.precise_chunks
    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size

    # Read BED/GFF and return genes and relevant flanking regions for calculating B
    blockstart, blockend =  \
        bedgffHandler(file_path) 

    # Initialize the array for B values (all initially set to 1.0)
    b_values = np.ones(chr_end - chr_start, dtype=np.float64)

    # Calculate cumulative conserved length in each chunk 
    lperchunk = calculateLPerChunk(chunk_size, blockstart, blockend, chr_start, chr_end)


    # # Iterate over chunks, calculating B for all neutral sites
    # for chunk_num in range(num_chunks): #Iterate through each chunk (Old loop)
    #     b_values = \
    #     process_single_chunk(chunk_num, chunk_size, blockstart, blockend, chr_start, chr_end, num_chunks, precise_chunks, lperchunk, b_values)
    for s, e in zip(blockstart, blockend): # Converts gene sites to NaN
        b_values[s - chr_start : e - chr_start + 1] = np.nan

    with ThreadPoolExecutor() as executor:
        results = [executor.submit(process_single_chunk, x, chunk_size, blockstart, blockend, chr_start, 
                                   chr_end, num_chunks, precise_chunks, lperchunk, b_values)
            for x in range(num_chunks)]
        
    print("Mean B (no nan sites):", b_values[~np.isnan(b_values)].mean())
    print("Total conserved length: ", int(sum(lperchunk)), round((sum(lperchunk)/(chr_end - chr_start))*100,2), "%" )

    if args.pop_change:
        return b_values#get_Bcur(b_values)
    else:
        return b_values
    
        # # Converts B at gene sites to 'nan'
    # for start, end in zip(blockstart, blockend):
    #     sites['B'][start - sites['pos'][0]:end - sites['pos'][0]] = np.nan 

    # nogene_sites = np.sort(sites[~np.isnan(sites['B'])], order=['B']) # Removes gene sites
    # print(nogene_sites)