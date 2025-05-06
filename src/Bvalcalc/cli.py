#!/usr/bin/env python3
import time
import os
import argparse
from bvalcalc.core.utils.parseArgs import parseGenomeArgs, parseRegionArgs, parseGeneArgs, parseSiteArgs
from bvalcalc.core.plotB import plotB
from bvalcalc.core.utils.generateParams import SPECIES, generateParams
import sys

def main():
    start_time = time.time()

    # check_generate_params_args() # Unique error message for --generate_params to print species names
    parser = argparse.ArgumentParser(description="Bcalc main function! :p")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--generate_params',metavar='SPECIES',nargs='?',const='template',default=None,choices=['human', 'drosophila', 'arabidopsis', 'mouse', 'pfalciparum', 'celegans', 'template'], help="Generate popgen params for a given species (human, drosophila, arabidopsis or mouse)")
    parser.add_argument('--dir', '-d', default='.', help="Directory in which to write the generated params file (default: current directory)")
    group.add_argument('--genome', '-w', action='store_true', help="Compute B values genome-wide for all sites considering all selected elements")
    group.add_argument('--region', '-r', action='store_true', help="Compute B values for a specific chromosomal region, considering genome-wide effects")
    group.add_argument('--gene', '-g', action='store_true', help="Compute B values for a neutral region adjacent to a single selected element")
    group.add_argument('--site', '-s',action='store_true', help="Compute B values for a single site from a selected element")
    known_args, remaining_args = parser.parse_known_args()

    if known_args.generate_params is not None: # if --generate_params
        print(f"Retrieving params from template...")
        generateParams(known_args.generate_params, known_args.dir)
        return

    print(f"= Calculating relative diversity (B) for all neutral sites across the genome. = = =")

    if known_args.genome: # Run genome Bcalc
        args = parseGenomeArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        from bvalcalc.core.genomeBcalc import genomeBcalc
        genomeBcalc(args)

    elif known_args.region: # Run region Bcalc
        args = parseRegionArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        from bvalcalc.core.regionBcalc import regionBcalc
        output_data, block_ranges = regionBcalc(args)
        if getattr(args, 'plot_output', True):
            plotB(b_values_input=output_data, caller="chromosome", output_path=args.plot_output, quiet=args.quiet, gene_ranges=block_ranges, neutral_only=args.neutral_only)

    elif known_args.gene: # Run gene Bcalc
        args = parseGeneArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        from bvalcalc.core.geneBcalc import geneBcalc
        output_data = geneBcalc(args) # Capture the output from geneBcalc
        if getattr(args, 'plot_output', False): # If the --plot_basic flag was provided, call plotB with geneBcalc's output.
            plotB(b_values_input=output_data, caller="gene", output_path=args.plot_output, quiet=args.quiet)

    elif known_args.site: # Run single site Bcalc
        args = parseSiteArgs(remaining_args)
        os.environ["BCALC_POP_PARAMS"] = args.pop_params  # Handle Params file
        from bvalcalc.core.siteBcalc import siteBcalc
        siteBcalc(args)

    print(f"= B value calculated in {time.time() - start_time:.2f} seconds. = = =")

if __name__ == "__main__":
    main()