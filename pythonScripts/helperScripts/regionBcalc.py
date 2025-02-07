from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.calculateB import calculateB
import numpy as np

def regionBcalc(args):    
    gene_size, flank_len = args.gene_size, args.flank_len

    b_values = calculateB(np.arange(1, flank_len, 1, dtype = int), gene_size)

    print(f"B for adjacent site: {b_values[0]}")
    print(f"Mean B for flanking region: {b_values.mean()}")
    print(f"All B: {b_values}")

    if args.pop_change:
        return get_Bcur(b_values)
    else:
        return b_values