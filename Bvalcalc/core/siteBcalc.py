from Bvalcalc.core.helpers.demography_helpers import get_Bcur

def siteBcalc(args):    
    gene_size, distance, quiet = args.gene_size, args.distance, args.quiet

    import Bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    dfe_helper.CONSTANT_DFE = args.constant_dfe # Update DFE if --constant_dfe
    from Bvalcalc.core.calculateB import calculateB_linear

# ## To delete
#     from Bvalcalc.core.calculateB import calculateB_hri


#     ## r = 1e-8
#     ## chunk_size = 20000
#     ## First we need to identify our linkage block of interest, which will be the area within a cM recombint length region, though when its shorter than 
#     ## linkage_block_max_rec_length = r * chunk_size # r * l
#     ## linkage_block_min_physical_length = 100000


#     ### THEN RUN
#     interfering_L = 10000
#     prior_B = 0.9

#     b_values = calculateB_hri(interfering_L, prior_B)
#     return
# ## To delete

    b_values = calculateB_linear(distance, gene_size)
    print(f"B for site {distance}bp away from {gene_size}bp region: {b_values}")

    if args.pop_change:
        if not quiet: print("Demographic change prior to B-calculation", b_values)
        b_values = get_Bcur(b_values)
        if not quiet: print("Demographic change applied to B-calculation", b_values)
    
    return