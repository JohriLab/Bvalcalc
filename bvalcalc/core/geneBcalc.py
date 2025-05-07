from bvalcalc.core.helpers.demography_helpers import get_Bcur
from bvalcalc.utils.bin_outputs import bin_outputs
import numpy as np
import os

def geneBcalc(args):    
    gene_size, flank_len, quiet = args.gene_size, args.flank_len, args.quiet
    
    import bvalcalc.utils.dfeHelper as dfeHelper
    dfeHelper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    from bvalcalc.core.calculateB import calculateB_linear # Import calculateB with updated gamma DFE if needed

    if not quiet: 
        print(f"====== P A R A M E T E R S =========================")
        print(f"Distribution of fitness effects (DFE): {flank_len}bp")
        print(f"Length of element under selection: {gene_size}bp")
        print(f"Length of flanking neutral region: {flank_len}bp")

    print(f"====== S T A R T I N G ===== C A L C ===============")
    b_values = calculateB_linear(np.arange(1, flank_len, 1, dtype = int), gene_size)
    print(f"====== F I N I S H E D ===== C A L C ===============")

    if not quiet:
        print(f"====== R E S U L T S ! =============================")

    if args.pop_change:
        if not quiet: print("B prior to demographic change", b_values)
        b_values = get_Bcur(b_values)
        if not quiet: print("B post B-calculation", b_values)
    output_data = np.column_stack((np.arange(1, flank_len, 1, dtype = int), b_values))

    if not quiet:
        print(f"B for adjacent site: {b_values[0]}")
        print(f"Mean B for flanking region: {b_values.mean()}")
        print(f"B at start and end of the neutral region: {b_values}")

    if args.out is not None: # Write to CSV
        # 1) decouple
        positions = output_data[:, 0].astype(int)
        bvals     = output_data[:, 1].astype(float)
        b_bvals, b_pos = bin_outputs(bvals, positions, args.out_binsize) # 2) bin them
        binned_output = np.column_stack((b_pos, b_bvals)) # 3) rebuild a two-column array for the binned output
        
        np.savetxt(args.out, # This might be "b_values.csv" or a custom path
            binned_output, delimiter=",", header="Distance,B", fmt=("%d", "%.6f"), comments="")
        print(f"Saved B values to: {os.path.abspath(args.out)}")
    else:
        if not args.quiet:
            print("No output CSV requested; skipping save.")
    
    return output_data