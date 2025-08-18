from Bvalcalc.core.helpers.demography_helpers import get_Bcur

def siteBcalc(args):    
    element_size, distance, quiet = args.element_size, args.distance, args.quiet

    import Bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    dfe_helper.CONSTANT_DFE = args.constant_dfe # Update DFE if --constant_dfe
    from Bvalcalc.core.calculateB import calculateB_linear

    b_values = calculateB_linear(distance, element_size)
    print(f"B for site {distance}bp away from {element_size}bp region: {b_values}")

    if args.pop_change:
        if not quiet: print("Demographic change prior to B-calculation", b_values)
        b_values = get_Bcur(b_values)
        if not quiet: print("Demographic change applied to B-calculation", b_values)
    
    return