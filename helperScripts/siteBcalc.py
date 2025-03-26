from helperScripts.RunBCalcScripts.demographyHelpers import get_Bcur
from helperScripts.calculateB import calculateB_linear

def siteBcalc(args):    
    gene_size, distance, silent = args.gene_size, args.distance, args.silent

    b_values = calculateB_linear(distance, gene_size)
    print(f"B for site {distance}bp away from {gene_size}bp region: {b_values}")

    if args.pop_change:
        if not silent: print("Demographic change prior to B-calculation", b_values)
        b_values = get_Bcur(b_values)
        if not silent: print("Demographic change applied to B-calculation", b_values)
    
    return