from helperScripts.RunBCalcScripts.process_single_chunk import process_single_chunk
from helperScripts.RunBCalcScripts.bedgffHandler import bedgffHandler
from helperScripts.RunBCalcScripts.calculateLPerChunk import calculateLPerChunk
from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.RunBCalcScripts.recmapHandler import recmapHandler, calcRLengths
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

#Main function
def genomeBcalc(args):    
    file_path, chr_start, chr_end, chunk_size, precise_chunks = args.file_path, args.chr_start, args.chr_end, args.chunk_size, args.precise_chunks

    print(f"Calculating relative diversity (B) for all neutral sites across the genome...")

    print(f"====== P A R A M E T E R S =========================")
    print(f"BED/GFF file for regions under selection: {file_path}")
    print(f"First position in chromosome: {chr_start}")
    print(f"Last position in chromosome: {chr_end}")
    print(f"Size of chunks to calculate B in per iteration: {chunk_size}bp")
    print(f"Number of adjacent chunks to calculate B precisely for: {precise_chunks}")


    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size

    # Read BED/GFF and return genes and relevant flanking regions for calculating B
    blockstart, blockend =  \
        bedgffHandler(file_path) 

        # Process recombination map if provided
    if args.rec_map:
        print(f"Using recombination map from {args.rec_map}")
        rec_rate_per_chunk = recmapHandler(args.rec_map, chr_start, chr_end, chunk_size)
    else:
        rec_rate_per_chunk = None

    blockRLengths = calcRLengths(blockstart, blockend, rec_rate_per_chunk, chr_end, chr_start, chunk_size)

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

    print(f"====== R E S U L T S == P E R == C H U N K =========")
    with ThreadPoolExecutor() as executor:
        results = [executor.submit(process_single_chunk, x, chunk_size, blockstart, blockend, chr_start, 
                                   chr_end, num_chunks, precise_chunks, lperchunk, b_values, rec_rate_per_chunk)
            for x in range(num_chunks)]
    
    print(f"====== R E S U L T S ====== S U M M A R Y ===========")
    print(f"Mean B of neutral sites across genome: {b_values[~np.isnan(b_values)].mean()}bp")
    print(f"Cumulative length of regions under selection: {int(sum(lperchunk))}bp ({round((sum(lperchunk)/(chr_end - chr_start))*100,2)}%)")


    if args.pop_change:
        return b_values#get_Bcur(b_values)
    else:
        return b_values
    
        # # Converts B at gene sites to 'nan'
    # for start, end in zip(blockstart, blockend):
    #     sites['B'][start - sites['pos'][0]:end - sites['pos'][0]] = np.nan 

    # nogene_sites = np.sort(sites[~np.isnan(sites['B'])], order=['B']) # Removes gene sites
    # print(nogene_sites)