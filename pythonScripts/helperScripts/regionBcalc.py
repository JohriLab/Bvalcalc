from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.calculateB import calculateB
import numpy as np
import os

def regionBcalc(args):    
    gene_size, flank_len = args.gene_size, args.flank_len
    print(f"Calculating relative diversity (B) for a neutral region adjacent to a single selected region (e.g. gene, exon, regulatory element)...")

    print(f"====== P A R A M E T E R S =========================")
    print(f"Distribution of fitness effects (DFE): {flank_len}bp")
    print(f"Length of region under selection: {gene_size}bp")
    print(f"Length of flanking neutral region: {flank_len}bp")

    b_values = calculateB(np.arange(1, flank_len, 1, dtype = int), gene_size, rdistance_to_element=None, rlength_of_element=None)

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
        print("No output CSV requested; skipping save.")

    if args.pop_change:
        return get_Bcur(b_values)
    else:
        return b_values