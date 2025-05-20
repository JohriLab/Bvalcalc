import sys
import argparse

def parseSiteArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for a single neutral site given a distance from a single selected region and prints to console.")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--gene_size', type=int, default=10000, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--distance', type=int, default=1, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--gamma_dfe', action='store_true', help="If set, gamma distribution parameters will be used to define DFE's discretized f0-f3 proportions")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    return parser.parse_args(argv)

def parseGeneArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for neutral sites flanking a single element under selection.")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--gene_size', type=int, default=10000, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--flank_len', type=int, default=40000, help="Length of flanking neutral region for which to calcuate recovery of B. [25000]")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--gamma_dfe', action='store_true', help="If set, gamma distribution parameters will be used to define DFE's discretized f0-f3 proportions")
    parser.add_argument('--plot_output', nargs='?', const='Bplot.png', default=None, 
                        help="Generate a B recovery slope output"
                            "Provide path to plot output.")
    parser.add_argument('--out', type=str, default=None,
                        help="Optional path to output CSV file.")
    parser.add_argument('--out_binsize', type=int, default=None, help="Size of bins to write average B in. By default B is saved per-base")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    raw = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv)

    if '--out' in raw and '--out_binsize' not in raw: # enforce: if they asked for --out, they must have also provided --out_binsize
        parser.error("argument --out_binsize is required when --out is specified")
    return args

def parseRegionArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    # parser.add_argument('--pop_params', type=int, required=True, help="Path to file providing popgen parameters specific to modelled population (empirical or simulated).")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--bedgff_path', type=str, required=True, help="Path to input BED or GFF3 file.")
    parser.add_argument('--chunk_size', type=int, default=20000, help="Size of chunks calculated simultaneously (bp). [100000]")
    parser.add_argument('--precise_chunks', type=int, default=3, help="Number of adjacent chunks to calculate B precisely.")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--prior_Bmap', type=str, default=None,
                        help="Optional input with per-site expected diversity, e.g. a B map calculated using Bvalcalc on different annotations! Usage: --prior_Bmap your_map.csv "
                             "These values will be multiplied by the caluculated B, i.e. a value of 0.9 at a given position will be returned as 0.9 * [newly calculated B]"
                             "Format should be the same as the B-map output: 'Chromosome,Position,Conserved,B'. "
                             "Note that the Conserved column is needed for accurate parsing but will not affect the analysis.")
    parser.add_argument('--rec_map', nargs='?', default=None,
                        help="Optional recombination (crossover) map input. Usage: --rec_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that recombination rates will be averaged within each chunk.")    
    parser.add_argument('--gc_map', nargs='?', default=None,
                        help="Optional gene conversion (non-crossover) map input. Usage: --gc_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that gene conversion rates will be averaged within each chunk.")    
    parser.add_argument('--gamma_dfe', action='store_true', help="If set, gamma distribution parameters will be used to define DFE's discretized f0-f3 proportions")
    parser.add_argument('--plot_output', nargs='?', const='genome_plot.png', default=None, 
                        help="Generate a basic plot using `Bvalcalc.py --genome` output"
                            "Provide path to plot output.")
    parser.add_argument('--neutral_only', action='store_true', help="If set, plot_output will only show neutral sites.")
    parser.add_argument('--out', type=str, default=None,
                        help="Required path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved. Note that by default it is per-base B, to output B averaged across"
                             "bins, use --out_bins [int]")
    parser.add_argument('--out_binsize', type=int, default=None, help="Size of bins to write average B in. By default B is saved per-base")
    parser.add_argument('--verbose', action='store_true', help="If set, will give per-chunk summaries")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    
    raw = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv)

    if '--out' in raw and '--out_binsize' not in raw: # enforce: if they asked for --out, they must have also provided --out_binsize
        parser.error("argument --out_binsize is required when --out is specified")

        
    return args

def parseGenomeArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    # parser.add_argument('--pop_params', type=int, required=True, help="Path to file providing popgen parameters specific to modelled population (empirical or simulated).")
    parser.add_argument('--pop_params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--bedgff_path', type=str, required=True, help="Path to input BED or GFF3 file.")
    parser.add_argument('--chr_sizes', type=str, default=None, help="Chromosome sizes file. Defaults to end of last gene in each chromosome if not provided.")
    parser.add_argument('--chunk_size', type=int, default=20000, help="Size of chunks calculated simultaneously (bp). [100000]")
    parser.add_argument('--precise_chunks', type=int, default=3, help="Number of adjacent chunks to calculate B precisely.")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--prior_Bmap', type=str, default=None,
                        help="Optional input with per-site expected diversity, e.g. a B map calculated using Bvalcalc on different annotations! Usage: --prior_Bmap your_map.csv "
                             "These values will be multiplied by the caluculated B, i.e. a value of 0.9 at a given position will be returned as 0.9 * [newly calculated B]"
                             "Format should be the same as the B-map output: 'Chromosome,Position,Conserved,B'. "
                             "Note that the Conserved column is needed for accurate parsing but will not affect the analysis.")
    parser.add_argument('--rec_map', type=str, default=None,
                        help="Optional recombination (crossover) map input. Usage: --rec_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that recombination rates will be averaged within each chunk.")    
    parser.add_argument('--gc_map', type=str, default=None,
                        help="Optional gene conversion (non-crossover) map input. Usage: --gc_map your.map, "
                             "Format should be a two column csv with the header: 'start,rate'. "
                             "Note that gene conversion rates will be averaged within each chunk.") 
    parser.add_argument('--gamma_dfe', action='store_true', help="If set, gamma distribution parameters will be used to define DFE's discretized f0-f3 proportions")   
    parser.add_argument('--neutral_only', action='store_true', help="If set, plot_output will only show neutral sites.")
    parser.add_argument('--out', type=str, default=None,
                        help="Required path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved.")
    parser.add_argument('--out_binsize', type=int, default=None, help="Size of bins to write average B in. By default B is saved per-base")
    parser.add_argument('--verbose', action='store_true', help="If set, will give per-chunk summaries")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")

    raw = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv)

    if '--out' in raw and '--out_binsize' not in raw: # enforce: if they asked for --out, they must have also provided --out_binsize
        parser.error("argument --out_binsize is required when --out is specified")
        
    return args

def parseVcfArgs(argv=None):
    parser = argparse.ArgumentParser(description="VCF utilities for getting B statistics and filtering using on a B-map.")
    parser.add_argument('--Bmap', type=str, required=True, help="Path to B map, e.g., allchroms_B.csv")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    return parser.parse_args(argv)

def parse_args():
    parser = argparse.ArgumentParser(description="Welcome to Bvalcalc! Please specify a mode to calculate B.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--generate_params', metavar='SPECIES', nargs='?', const='template', default=None,
                       choices=['selfing', 'human', 'drosophila', 'arabidopsis', 'mouse', 'pfalciparum', 'celegans', 'template'],
                       help="Save population parameters from a species template")
    parser.add_argument('--dir', '-d', default='.', help="Directory to write the generated params file (default: current directory)")
    group.add_argument('--site', '-s', action='store_true', help="Calculate B values for a single site from a selected element")
    group.add_argument('--gene', '-g', action='store_true', help="Calculate B values for a neutral region adjacent to a single selected element")
    group.add_argument('--region', '-r', type=str, default=None, help="Calculate B values for a specific chromosomal region, considering genome-wide effects. Provide region as following argument [CHR,START,END].")
    group.add_argument('--genome', '-w', action='store_true', help="Calculate B values genome-wide for all sites considering all selected elements")
    group.add_argument('--vcf', '-v', type=str, default=None, help="VCF utilities for getting B statistics and filtering using a B-map. Provide VCF input as following argument.")
    return parser
