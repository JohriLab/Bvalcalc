#!/usr/bin/env python3
import os
import sys
import time
import argparse
from bvalcalc.utils.parseArgs import parse_args, parseGenomeArgs, parseRegionArgs, parseGeneArgs, parseSiteArgs, parseBmapArgs
from bvalcalc.core.plotB import plotB
from bvalcalc.utils.generateParams import SPECIES, generateParams, check_generate_params_args
from bvalcalc.core.positionsBstats import positionsBstats
from bvalcalc.core.plotChromB import plotChromB

def main():
    start_time = time.time()

    check_generate_params_args() # Unique error message for --generate_params to print species names
    parser = parse_args()
    known_args, remaining_args = parser.parse_known_args()

    if known_args.generate_params is not None: # if --generate_params
        print(f"Retrieving params from template...")
        generateParams(known_args.generate_params, known_args.dir)
        return
    
    if known_args.Bmap is not None: # if --generate_params
        args = parseBmapArgs(remaining_args)
        flat_b, flat_chrom = positionsBstats(args, known_args.Bmap)
        if args.plot_distribution:
            plotChromB(flat_b, flat_chrom, args.plot_distribution, args.quiet)
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
        output_data, block_ranges = regionBcalc(args, known_args.region)
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