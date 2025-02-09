from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.calculateB import calculateB
import numpy as np

def regionBcalc(args):    
    gene_size, flank_len = args.gene_size, args.flank_len
    print(f"Calculating relative diversity (B) for a neutral region adjacent to a single selected region (e.g. gene, exon, regulatory element)...")

    print(f"====== P A R A M E T E R S =========================")
    print(f"Distribution of fitness effects (DFE): {flank_len}bp")
    print(f"Length of region under selection: {gene_size}bp")
    print(f"Length of flanking neutral region: {flank_len}bp")

    b_values = calculateB(np.arange(1, flank_len, 1, dtype = int), gene_size)

    print(f"====== R E S U L T S ! =============================")
    print(f"B for adjacent site: {b_values[0]}")
    print(f"Mean B for flanking region: {b_values.mean()}")
    print(f"B at start and end of the neutral region: {b_values}")

    if args.pop_change:
        return get_Bcur(b_values)
    else:
        return b_values