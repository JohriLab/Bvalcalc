from core.chromBcalc import chromBcalc
from core.utils.bedgffHandler import bedgffHandler
import numpy as np
import sys

def regionBcalc(args):    

    allblockstart, allblockend, allblockchrom,  = bedgffHandler(args.bedgff_path) # Read BED/GFF, return start and end of conserved elements

    calc_chrom, start, end = parse_region(args.calc_region)

    # # First, require calc_region
    # # Then, mask for the chromosome the calc_region is in.
    # # Then, pass that into chromBcalc() as calc_start and calc_end
    mask = allblockchrom == calc_chrom
    blockstart = allblockstart[mask]
    blockend = allblockend[mask]
    chromosome = calc_chrom
    output_data, block_ranges = chromBcalc(args, blockstart, blockend, chromosome)

    # print("gaten", unique_chromosomes[i])
    sys.exit()


    # unique_chromosomes = np.unique(allblockchrom) # Move BED/GFF handler here
    # print(unique_chromosomes) ## Now, loop over each chromosome and save B output
    # for i in np.arange(0,len(unique_chromosomes)):
    #     mask = allblockchrom == unique_chromosomes[i]
    #     blockstart = allblockstart[mask]
    #     blockend = allblockend[mask]
    #     chromosome = unique_chromosomes[i]
    #     output_data, block_ranges = chromBcalc(args, blockstart, blockend, chromosome)
    
    # # Don't need to capture output, simply save in chromBcalc or not. 
    # # Then, if you want to plot a region, must use --calc_region chr:start-end (instead of calc_start,calc_end)
    # # Pass up plot data as needed

    return  output_data, block_ranges

def parse_region(region_str):
    try:
        chrom_part, pos_part = region_str.split(":")
        start_str, end_str = pos_part.split("-")
        chrom = chrom_part
        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))
        return chrom, start, end
    except ValueError:
        raise ValueError(f"Region format invalid: '{region_str}' should be like 'chr1:12345-67890'")