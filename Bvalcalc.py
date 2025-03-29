#!/usr/bin/env python3
import sys
import time
import numpy as np
import os
import argparse
from core.utils.parseArgs import parseGenomeArgs, parseRegionArgs, parseSiteArgs
from core.plotB import plotB

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Bcalc main function! :p")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--genome', action='store_true', help="Compute B values genome-wide")
    group.add_argument('--region', action='store_true', help="Compute B values for a region")
    group.add_argument('--site', action='store_true', help="Compute B values for a single site")
    known_args, remaining_args = parser.parse_known_args()

    if known_args.genome: # Run genome Bcalc
        args = parseGenomeArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        from core.genomeBcalc import genomeBcalc
        output_data, block_ranges = genomeBcalc(args)
        if getattr(args, 'plot_output', True):
            plotB(b_values_input=output_data, caller="genome", output_path=args.plot_output, silent=args.silent, gene_ranges=block_ranges, neutral_only=args.neutral_only)

    elif known_args.region: # Run region Bcalc
        args = parseRegionArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        from core.regionBcalc import regionBcalc
        output_data = regionBcalc(args) # Capture the output from regionBcalc
        if getattr(args, 'plot_output', False): # If the --plot_basic flag was provided, call plotB with regionBcalc's output.
            plotB(b_values_input=output_data, caller="region", output_path=args.plot_output, silent=args.silent)

    elif known_args.site: # Run single site Bcalc
        args = parseSiteArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        from core.siteBcalc import siteBcalc
        siteBcalc(args)

    print(f"= B value calculated in {time.time() - start_time:.2f} seconds. = = =")

if __name__ == "__main__":
    main()