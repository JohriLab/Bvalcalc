import sys
import argparse

def parseSiteArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for a single neutral site given a distance from a single selected region and prints to console.")
    parser.add_argument('--params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--element_size', type=int, default=10000, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--distance', type=int, default=1, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--gamma_dfe', action='store_true', help="If set, gamma distribution parameters will be used to define DFE's discretized f0-f3 proportions")
    parser.add_argument('--constant_dfe', action='store_true', help="If set, the constant `s` and `neu_prop` parameters in the params file will be used for all conserved regions instead of discretized f0-f3 proportions")   
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    return parser.parse_args(argv)

def parseGeneArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for neutral sites flanking a single element under selection.")
    parser.add_argument('--params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--element_size', type=int, default=10000, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--flank_len', type=int, default=40000, help="Length of flanking neutral region for which to calcuate recovery of B. [25000]")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--gamma_dfe', action='store_true', help="If set, gamma distribution parameters will be used to define DFE's discretized f0-f3 proportions")
    parser.add_argument('--constant_dfe', action='store_true', help="If set, the constant `s` and `neu_prop` parameters in the params file will be used for all conserved regions instead of discretized f0-f3 proportions")   
    parser.add_argument('--plot', nargs='?', const='Bplot.png', default=None, 
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
    # parser.add_argument('--params', type=int, required=True, help="Path to file providing popgen parameters specific to modelled population (empirical or simulated).")
    parser.add_argument('--params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--bedgff', type=str, required=True, help="Path to input BED or GFF3 file.")
    parser.add_argument('--chunk_size', type=int, default=None, help="Size of chunks calculated simultaneously (bp). Default: 20000, auto-adjusted for large datasets.")
    parser.add_argument('--chunk_batch_size', type=int, default=250, help="Number of chunks to process before flushing memory. [250]")
    parser.add_argument('--precise_chunks', type=int, default=3, help="Number of adjacent chunks to calculate B precisely.")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--prior_Bmap', type=str, default=None,
                        help="Optional input with per-site expected diversity, e.g. a B map calculated using Bvalcalc on different annotations! Usage: --prior_Bmap your_map.csv "
                             "These values will be multiplied by the caluculated B, i.e. a value of 0.9 at a given position will be returned as 0.9 * [newly calculated B]"
                             "Format should be the same as the B-map output: 'Chromosome,Start,Conserved,B'. "
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
    parser.add_argument('--constant_dfe', action='store_true', help="If set, the constant `s` and `neu_prop` parameters in the params file will be used for all conserved regions instead of discretized f0-f3 proportions")   
    parser.add_argument('--plot', nargs='?', const='genome_plot.png', default=None, 
                        help="Generate a basic plot using `Bvalcalc --genome` output"
                            "Provide path to plot output.")
    parser.add_argument('--hri', action='store_true', help="If set, will enable post-hoc calculation of B under HRI (B'; Becher and Charlesworth 2025), for low recombination regions")   
    parser.add_argument('--neutral_only', action='store_true', help="If set, plot will only show neutral sites.")
    parser.add_argument('--out', type=str, default=None,
                        help="Required path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved. Note that by default it is per-base B, to output B averaged across"
                             "bins, use --out_bins [int]")
    parser.add_argument('--out_binsize', type=int, default=None, help="Size of bins to write average B in. By default B is saved per-base")
    parser.add_argument('--verbose', action='store_true', help="If set, will give per-chunk summaries")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    parser.add_argument('--no_header', action='store_true', help="If set, skip writing comment headers to output file (version, command, etc.)")
    
    raw = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv)

    if '--out' in raw and '--out_binsize' not in raw: # enforce: if they asked for --out, they must have also provided --out_binsize
        parser.error("argument --out_binsize is required when --out is specified")

        
    return args

def parseGenomeArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    # parser.add_argument('--params', type=int, required=True, help="Path to file providing popgen parameters specific to modelled population (empirical or simulated).")
    parser.add_argument('--params', type=str, required=True, help="Path to Python file with population genetic parameters, e.g., ExampleParams.py")
    parser.add_argument('--bedgff', type=str, required=True, help="Path to input BED or GFF3 file.")
    parser.add_argument('--chr_sizes', type=str, default=None, help="Chromosome sizes file. Defaults to end of last gene in each chromosome if not provided.")
    parser.add_argument('--chunk_size', type=int, default=None, help="Size of chunks calculated simultaneously (bp). Default: 20000, auto-adjusted for large datasets.")
    parser.add_argument('--chunk_batch_size', type=int, default=250, help="Number of chunks to process before flushing memory. [250]")
    parser.add_argument('--precise_chunks', type=int, default=3, help="Number of adjacent chunks to calculate B precisely.")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--prior_Bmap', type=str, default=None,
                        help="Optional input with per-site expected diversity, e.g. a B map calculated using Bvalcalc on different annotations! Usage: --prior_Bmap your_map.csv "
                             "These values will be multiplied by the caluculated B, i.e. a value of 0.9 at a given position will be returned as 0.9 * [newly calculated B]"
                             "Format should be the same as the B-map output: 'Chromosome,Start,Conserved,B'. "
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
    parser.add_argument('--constant_dfe', action='store_true', help="If set, the constant `s` and `neu_prop` parameters in the params file will be used for all conserved regions instead of discretized f0-f3 proportions")   
    parser.add_argument('--hri', action='store_true', help="If set, will enable post-hoc calculation of B under HRI (B'; Becher and Charlesworth 2025), for low recombination regions")   
    parser.add_argument('--neutral_only', action='store_true', help="If set, plot will only show neutral sites.")
    parser.add_argument('--out', type=str, default=None,
                        help="Required path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved.")
    parser.add_argument('--out_binsize', type=int, default=None, help="Size of bins to write average B in. By default B is saved per-base")
    parser.add_argument('--verbose', action='store_true', help="If set, will give per-chunk summaries")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    parser.add_argument('--no_header', action='store_true', help="If set, skip writing comment headers to output file (version, command, etc.)")

    raw = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv)

    if '--out' in raw and '--out_binsize' not in raw: # enforce: if they asked for --out, they must have also provided --out_binsize
        parser.error("argument --out_binsize is required when --out is specified")
        
    return args

def parseBmapArgs(argv=None):
    parser = argparse.ArgumentParser(description="B-map utilities for getting B statistics for sites in a VCF/txt file with specific positions.")
    parser.add_argument('--positions', type=str, required=True, help="VCF or  input as following argument")
    parser.add_argument('--plot_distribution', nargs='?', const='B_distribution.png', default=None, help="Output path for a plot of the distribution of B across each chromosome.")   
    parser.add_argument('--out', type=str, default=None,
                        help="Path to save per-site B for variant sites in the VCF/txt file. Results are not saved if --out is not specified.")
    parser.add_argument('--out_minimum', type=float, default=None, help="If set, only the sites ABOVE the given threshold of B will be returned, i.e. B > [threshold].")   
    parser.add_argument('--out_maximum', type=float, default=None, help="If set, only the sites BELOW the given threshold of B will be returned, i.e. B < [threshold].")   
    parser.add_argument('--bcftools_format', action='store_true', help="If set, removes the B value column and reformats the output for filter with bcftools view.")
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    return parser.parse_args(argv)

def parse_args(version):
    parser = argparse.ArgumentParser(
        prog="Bvalcalc",
        description=f"Welcome to Bvalcalc v{version}! Please specify a mode to calculate B."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--generate_params', metavar='SPECIES', nargs='?', const='template', default=None,
                       help="Save population parameters from a species template")
    group.add_argument('--site', '-s', action='store_true', help="Calculate B values for a single site from a selected element")
    group.add_argument('--gene', '-g', action='store_true', help="Calculate B values for a neutral region adjacent to a single selected element")
    group.add_argument('--region', '-r', type=str, help="Calculate B values for a specific chromosomal region, considering genome-wide effects. Provide region as [CHR,START,END].")
    group.add_argument('--genome', '-w', action='store_true', help="Calculate B values genome-wide for all sites considering all selected elements")
    group.add_argument('--Bmap', '-b', type=str, help="B-map lookup for sites in a VCF/txt file")
    group.add_argument('--download_sample_data', action='store_true', help="Download sample data files to current directory (or use --dir to specify)")
    
    parser.add_argument('--dir', '-d', default='.', help="Directory to write the generated params file or download sample data (default: current directory)")

    # If no args provided, print help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser
