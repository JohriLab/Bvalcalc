import argparse

def parseGenomeArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for all neutral sites across given chromosome.")
    parser.add_argument('--chr_start', type=int, required=True, help="Start of chromosome range.")
    parser.add_argument('--chr_end', type=int, required=True, help="End of chromosome range.")
    parser.add_argument('--chunk_size', type=int, default=25000, help="Size of chunks calculated simultaneously (bp). [100000]")
    parser.add_argument('--flank_len', type=int, default=100000, help="Length of region adjacent to conserved element (bp). [20000]")
    parser.add_argument('--file_path', type=str, required=True, help="Path to input BED or GFF3 file.")
    parser.add_argument('--precise_chunks', type=int, default=3, help="Number of adjacent chunks to calculate B precisely.")
    parser.add_argument('--pop_change', action='store_false', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    return parser.parse_args(argv)

def parseRegionArgs(argv=None):
    parser = argparse.ArgumentParser(description="Calculates B for neutral sites flanking a single region under selection.")
    parser.add_argument('--gene_size', type=int, default=10000, help="Length of single region (e.g. gene) under selection. [5000]")
    parser.add_argument('--flank_len', type=int, default=40000, help="Length of flanking neutral region for which to calcuate recovery of B. [25000]")
    parser.add_argument('--pop_change', action='store_true', help="If set, B will reflect the current B after a step change in population size, rather than ancestral B.")
    parser.add_argument('--plotBasic', action='store_true', help="Generate a basic plot using `Bvalcalc.py --region` output")
    parser.add_argument('--out', nargs='?', const='../../bin/b_values.csv', default=None,
                        help="Optional path to output CSV file. If --out is specified but no file name is given, "
                             "'b_values.csv' will be used in the current directory. If --out is not specified, "
                             "no CSV will be saved.")
    return parser.parse_args(argv)