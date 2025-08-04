from Bvalcalc.core.helpers.demography_helpers import get_Bcur

def siteBcalc(args):    
    gene_size, distance, quiet = args.gene_size, args.distance, args.quiet

    import Bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    dfe_helper.CONSTANT_DFE = args.constant_dfe # Update DFE if --constant_dfe
    from Bvalcalc.core.calculateB import calculateB_linear
    from Bvalcalc.core.calculateB import calculateB_hri
    # f1,f2,u,t1,t2,t3,N0,h
    interfering_L = 10000

    b_values = calculateB_hri(interfering_L)
    # b_values = calculateB_linear(distance, gene_size)
    print(f"B for site {distance}bp away from {gene_size}bp region: {b_values}")

    if args.pop_change:
        if not quiet: print("Demographic change prior to B-calculation", b_values)
        b_values = get_Bcur(b_values)
        if not quiet: print("Demographic change applied to B-calculation", b_values)
    
    return