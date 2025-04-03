from core.helpers.demography_helpers import get_Bcur
from core.calculateB import calculateB_linear
import numpy as np
import os

def regionBcalc(args):    
    gene_size, flank_len, quiet = args.gene_size, args.flank_len, args.quiet
    print(f"= Calculating relative diversity (B) for a neutral region adjacent to a single selected element = = =")
    if not quiet: 
        print(f"====== P A R A M E T E R S =========================")
        print(f"Distribution of fitness effects (DFE): {flank_len}bp")
        print(f"Length of region under selection: {gene_size}bp")
        print(f"Length of flanking neutral region: {flank_len}bp")

    print(f"====== S T A R T I N G ===== C A L C ===============")
    b_values = calculateB_linear(np.arange(1, flank_len, 1, dtype = int), gene_size)
    print(f"====== F I N I S H E D ===== C A L C ===============")

    if not quiet:
        print(f"====== R E S U L T S ! =============================")
        print(f"B for adjacent site: {b_values[0]}")
        print(f"Mean B for flanking region: {b_values.mean()}")
        print(f"B at start and end of the neutral region: {b_values}")

    if args.pop_change:
        if not quiet: print("Demographic change prior to B-calculation", b_values)
        b_values = get_Bcur(b_values)
        if not quiet: print("Demographic change applied to B-calculation", b_values)
    output_data = np.column_stack((np.arange(1, flank_len, 1, dtype = int), b_values))

    if args.out is not None: # Write to CSV
        np.savetxt(args.out, # This might be "b_values.csv" or a custom path
            output_data, delimiter=",", header="Distance,B", fmt=("%d", "%.6f"), comments="")
        print(f"Saved B values to: {os.path.abspath(args.out)}")
    else:
        if not args.quiet:
            print("No output CSV requested; skipping save.")
    
    return output_data