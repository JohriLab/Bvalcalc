from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.calculateB import calculateB_linear
import numpy as np
import os

def regionBcalc(args):    
    gene_size, flank_len, silent = args.gene_size, args.flank_len, args.silent
    print(f"= Calculating relative diversity (B) for a neutral region adjacent to a single selected element = = =")
    if not silent: 
        print(f"====== P A R A M E T E R S =========================")
        print(f"Distribution of fitness effects (DFE): {flank_len}bp")
        print(f"Length of region under selection: {gene_size}bp")
        print(f"Length of flanking neutral region: {flank_len}bp")

    print(f"====== S T A R T I N G ===== C A L C ===============")
    b_values = calculateB_linear(np.arange(1, flank_len, 1, dtype = int), gene_size)

    if not silent:
        print(f"====== R E S U L T S ! =============================")
        print(f"B for adjacent site: {b_values[0]}")
        print(f"Mean B for flanking region: {b_values.mean()}")
        print(f"B at start and end of the neutral region: {b_values}")

    if args.out is not None:
        csv_file = args.out  # This might be "b_values.csv" or a custom path

        # Combine positions and b_values into two columns
        output_data = np.column_stack((np.arange(1, flank_len, 1, dtype = int), b_values))

        # Write to CSV
        np.savetxt(
            csv_file, 
            output_data,
            delimiter=",", 
            header="Distance,B", 
            fmt=("%d", "%.6f"),  # first column = integer, second = float w/ 6 decimals
            comments=""
)
        print(f"Saved B values to: {os.path.abspath(csv_file)}")
    else:
        if not silent:
            print("No output CSV requested; skipping save.")

    if args.pop_change:
        return get_Bcur(b_values)
    else:
        return b_values