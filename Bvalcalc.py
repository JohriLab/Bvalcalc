#!/usr/bin/env python3
import sys
import time
import numpy as np
import os
import argparse
from core.utils.parseArgs import parseGenomeArgs, parseRegionArgs, parseGeneArgs, parseSiteArgs
from core.plotB import plotB

def main():
    start_time = time.time()

    print(f"= Calculating relative diversity (B) for all neutral sites across the genome. = = =")

    parser = argparse.ArgumentParser(description="Bcalc main function! :p")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--genome', action='store_true', help="Compute B values genome-wide for all sites considering all selected elements")
    group.add_argument('--region', action='store_true', help="Compute B values for a specific chromosomal region, considering genome-wide effects")
    group.add_argument('--gene', action='store_true', help="Compute B values for a neutral region adjacent to a single selected element")
    group.add_argument('--site', action='store_true', help="Compute B values for a single site from a selected element")
    known_args, remaining_args = parser.parse_known_args()

    if known_args.genome: # Run genome Bcalc
        args = parseGenomeArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        import core.utils.dfeHelper as dfeHelper
        dfeHelper.GAMMA_DFE = args.gamma_dfe
        from core.genomeBcalc import genomeBcalc
        genomeBcalc(args)

    elif known_args.region: # Run region Bcalc
        args = parseRegionArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        import core.utils.dfeHelper as dfeHelper
        dfeHelper.GAMMA_DFE = args.gamma_dfe
        from core.regionBcalc import regionBcalc
        output_data, block_ranges = regionBcalc(args)
        if getattr(args, 'plot_output', True):
            plotB(b_values_input=output_data, caller="chromosome", output_path=args.plot_output, quiet=args.quiet, gene_ranges=block_ranges, neutral_only=args.neutral_only)

    elif known_args.gene: # Run gene Bcalc
        args = parseGeneArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        import core.utils.dfeHelper as dfeHelper
        dfeHelper.GAMMA_DFE = args.gamma_dfe
        from core.geneBcalc import geneBcalc
        output_data = geneBcalc(args) # Capture the output from geneBcalc
        if getattr(args, 'plot_output', False): # If the --plot_basic flag was provided, call plotB with geneBcalc's output.
            plotB(b_values_input=output_data, caller="gene", output_path=args.plot_output, quiet=args.quiet)

    elif known_args.site: # Run single site Bcalc
        args = parseSiteArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        import core.utils.dfeHelper as dfeHelper
        dfeHelper.GAMMA_DFE = args.gamma_dfe
        from core.siteBcalc import siteBcalc
        siteBcalc(args)

    print(f"= B value calculated in {time.time() - start_time:.2f} seconds. = = =")

if __name__ == "__main__":
    main()