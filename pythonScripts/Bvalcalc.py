#!/usr/bin/env python3

# ./Bcalc_multiElements.py --calc_start 1 --calc_end 25254535 --file_path ../exampleData/dmel6_2R_genes.csv

#calculate_B_analytically_Eq3_mine_demography_multiple_elements
#This script is to get B values across a genomic element with multiple functional elements.
#Currently the output is in terms of an average of a sliding window
import sys
from helperScripts.parseArgs import parseGenomeArgs, parseRegionArgs
from helperScripts.genomeBcalc import genomeBcalc
from helperScripts.regionBcalc import regionBcalc
from helperScripts.plotBasic import plotBasic
import time
import numpy as np

import argparse

#CLI handling
def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Bcalc main! :p")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--genome', action='store_true', help="Compute B values genome-wide")
    group.add_argument('--region', action='store_true', help="Compute B values for a region")

    known_args, remaining_args = parser.parse_known_args()

    if known_args.genome:
        args = parseGenomeArgs(remaining_args)
        genome_output = genomeBcalc(args)
        if getattr(args, 'plotBasic', False):
            plotBasic(b_values_input=genome_output, caller="genome", silent=args.silent, output_file='../../bin/plotGenome.png')
    elif known_args.region:
        args = parseRegionArgs(remaining_args)
        # Capture the output from regionBcalc
        region_output = regionBcalc(args)
        # If the --plot_basic flag was provided, call plotBasic with regionBcalc's output.
        if getattr(args, 'plotBasic', False):
                plotBasic(b_values_input=region_output, caller="region", silent=args.silent, output_file='../../bin/plotRegion.png')

    print(f"= B value calculated in {time.time() - start_time:.2f} seconds. = = =")

if __name__ == "__main__":
    main()

sys.exit()

    # plotBasic(args)
    # siteBcalc(args)