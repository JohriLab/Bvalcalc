#!/usr/bin/env python3
import os
import sys
import time
import argparse
from Bvalcalc.utils.parseArgs import parse_args, parseGenomeArgs, parseRegionArgs, parseGeneArgs, parseSiteArgs, parseBmapArgs
from Bvalcalc.core.plotB import plotB
from Bvalcalc.core.deprecated.plotB_figures import plotB_figures
from Bvalcalc.core.deprecated.plotB_figures_200kb import plotB_figures_200kb
from Bvalcalc.utils.generateParams import SPECIES, generateParams, check_generate_params_args
from Bvalcalc.core.positionsBstats import positionsBstats
from Bvalcalc.core.plotChromB import plotChromB
from Bvalcalc.utils.sample_data import download_sample_data

__version__ = "1.3.0"

def main():
    start_time = time.time()

    check_generate_params_args() # Unique error message for --generate_params to print species names
    parser = parse_args(__version__)
    known_args, remaining_args = parser.parse_known_args()

    if known_args.generate_params is not None: # if --generate_params
        print(f"Retrieving params from template...")
        generateParams(known_args.generate_params, known_args.dir)
        return
    
    if known_args.download_sample_data: # if --download_sample_data
        success = download_sample_data(force=False, quiet=False, target_dir=known_args.dir)
        sys.exit(0 if success else 1)
    
    if known_args.Bmap is not None: # if --Bmap
        args = parseBmapArgs(remaining_args)
        flat_b, flat_chrom = positionsBstats(args, known_args.Bmap)
        if args.plot_distribution:
            plotChromB(flat_b, flat_chrom, args.plot_distribution, args.quiet)
        return

    print(f"= Calculating relative diversity (B) for all neutral sites across the genome. = = =")

    if known_args.genome: # Run genome Bcalc
        args = parseGenomeArgs(remaining_args)
        os.environ["BCALC_params"] = args.params  # Save params to global
        from Bvalcalc.core.genomeBcalc import genomeBcalc
        genomeBcalc(args)

    elif known_args.region: # Run region Bcalc
        args = parseRegionArgs(remaining_args)
        os.environ["BCALC_params"] = args.params  # Save params to global
        from Bvalcalc.core.regionBcalc import regionBcalc
        output_data, block_ranges, rec_rate_per_chunk_in_region, chunk_size = regionBcalc(args, known_args.region)
        if getattr(args, 'plot', True):
            plotB(b_values_input=output_data, caller="chromosome", output_path=args.plot, quiet=args.quiet, gene_ranges=block_ranges, neutral_only=args.neutral_only, rec_rates=rec_rate_per_chunk_in_region, chunk_size=chunk_size)

    elif known_args.gene: # Run gene Bcalc
        args = parseGeneArgs(remaining_args)
        os.environ["BCALC_params"] = args.params  # Save params to global
        from Bvalcalc.core.geneBcalc import geneBcalc
        output_data = geneBcalc(args) # Capture the output from geneBcalc
        if getattr(args, 'plot', False): # If the --plot flag was provided, call plotB with geneBcalc's output.
            plotB(b_values_input=output_data, caller="gene", output_path=args.plot, quiet=args.quiet)

    elif known_args.site: # Run single site Bcalc
        args = parseSiteArgs(remaining_args)
        os.environ["BCALC_params"] = args.params  # Save params to global
        from Bvalcalc.core.siteBcalc import siteBcalc
        siteBcalc(args)

    print(f"= B value calculated in {time.time() - start_time:.2f} seconds. = = =")

if __name__ == "__main__":
    main()