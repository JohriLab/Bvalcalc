import argparse

def parseGenomeArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    # parser.add_argument('--pop_params', type=int, required=True, help="Path to file providing popgen parameters specific to modelled population (empirical or simulated).")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--bedgff_path', type=str, required=True, help="Path to input BED or GFF3 file.")
    parser.add_argument('--chr_end', type=int, default=None, help="End of chromosome position. Defaults to calc_end if not provided.")
    parser.add_argument('--calc_start', type=int, default=None, help="Start of region to calculate B [default: chr_start]") # See if statment below
    parser.add_argument('--calc_end', type=int, default=None, help="End of region to calculate B [default: chr_end]") # See if statment below
    parser.add_argument('--chunk_size', type=int, default=20000, help="Size of chunks calculated simultaneously (bp). [100000]")
    parser.add_argument('--precise_chunks', type=int, default=3, help="Number of adjacent chunks to calculate B precisely.")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--rec_map', nargs='?', default=None,
                        help="Optional recombination (crossover) map input. Usage: --rec_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that recombination rates will be averaged within each chunk.")    
    parser.add_argument('--gc_map', nargs='?', default=None,
                        help="Optional gene conversion (non-crossover) map input. Usage: --gc_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that gene conversion rates will be averaged within each chunk.")    
    parser.add_argument('--plot_output', nargs='?', const='Bplot.png', default=None, 
                        help="Generate a basic plot using `Bvalcalc.py --genome` output"
                            "Provide path to plot output.")
    parser.add_argument('--neutral_only', action='store_true', help="If set, plot_output will only show neutral sites.")
    parser.add_argument('--out', nargs='?', default=None,
                        help="Required path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved.")
    parser.add_argument('--verbose', action='store_true', help="If set, will give per-chunk summaries")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    
    args = parser.parse_args(argv)
        
    return args

def parseRegionArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    # parser.add_argument('--pop_params', type=int, required=True, help="Path to file providing popgen parameters specific to modelled population (empirical or simulated).")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--bedgff_path', type=str, required=True, help="Path to input BED or GFF3 file.")
    parser.add_argument('--chr_end', type=int, default=None, help="End of chromosome position. Defaults to calc_end if not provided.")
    parser.add_argument('--calc_start', type=int, default=None, help="Start of region to calculate B [default: chr_start]") # See if statment below
    parser.add_argument('--calc_end', type=int, default=None, help="End of region to calculate B [default: chr_end]") # See if statment below
    parser.add_argument('--chunk_size', type=int, default=20000, help="Size of chunks calculated simultaneously (bp). [100000]")
    parser.add_argument('--precise_chunks', type=int, default=3, help="Number of adjacent chunks to calculate B precisely.")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--rec_map', nargs='?', default=None,
                        help="Optional recombination (crossover) map input. Usage: --rec_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that recombination rates will be averaged within each chunk.")    
    parser.add_argument('--gc_map', nargs='?', default=None,
                        help="Optional gene conversion (non-crossover) map input. Usage: --gc_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that gene conversion rates will be averaged within each chunk.")    
    parser.add_argument('--plot_output', nargs='?', const='genome_plot.png', default=None, 
                        help="Generate a basic plot using `Bvalcalc.py --genome` output"
                            "Provide path to plot output.")
    parser.add_argument('--neutral_only', action='store_true', help="If set, plot_output will only show neutral sites.")
    parser.add_argument('--out', nargs='?', default=None,
                        help="Required path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved.")
    parser.add_argument('--verbose', action='store_true', help="If set, will give per-chunk summaries")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    
    args = parser.parse_args(argv)
        
    return args

def parseGeneArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for neutral sites flanking a single element under selection.")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--gene_size', type=int, default=10000, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--flank_len', type=int, default=40000, help="Length of flanking neutral region for which to calcuate recovery of B. [25000]")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--plot_output', nargs='?', const='Bplot.png', default=None, 
                        help="Generate a basic plot using `Bvalcalc.py --genome` output"
                            "Provide path to plot output.")
    parser.add_argument('--out', nargs='?', default=None,
                        help="Optional path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved.")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    return parser.parse_args(argv)

def parseSiteArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for a single neutral site given a distance from a single selected region and prints to console.")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--gene_size', type=int, default=10000, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--distance', type=int, default=1, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    return parser.parse_args(argv)
