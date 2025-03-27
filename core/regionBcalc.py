from core.RunBCalcScripts.demographyHelpers import get_Bcur
from core.calculateB import calculateB_linear
import numpy as np

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
    print(f"====== F I N I S H E D ===== C A L C ===============")

    if not silent:
        print(f"====== R E S U L T S ! =============================")
        print(f"B for adjacent site: {b_values[0]}")
        print(f"Mean B for flanking region: {b_values.mean()}")
        print(f"B at start and end of the neutral region: {b_values}")

    if args.pop_change:
        if not silent: print("Demographic change prior to B-calculation", b_values)
        b_values = get_Bcur(b_values)
        if not silent: print("Demographic change applied to B-calculation", b_values)
    output_data = np.column_stack((np.arange(1, flank_len, 1, dtype = int), b_values))
    
    return output_data