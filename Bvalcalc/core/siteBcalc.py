from Bvalcalc.core.helpers.demography_helpers import get_Bcur

def siteBcalc(args):    
    gene_size, distance, quiet = args.gene_size, args.distance, args.quiet

    import Bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    dfe_helper.CONSTANT_DFE = args.constant_dfe # Update DFE if --constant_dfe
    from Bvalcalc.core.calculateB import calculateB_linear

# ## To delete
    from Bvalcalc.core.calculateB import calculateB_hri


    ## Get normal B

    ## r = 1e-8
    ## chunk_size = 20000
    ## First we need to identify our linkage block of interest, which will be the area within a cM recombint length region, with a minimum U

    ##      0. Make a normal B-map.

    ##      1. Find consecutive chunks with local rec_rate below 0.1 * r with int_L > 0 to have interference calculation applied within them.
    ## 
    ##      2. Re-calculate B from everything EXCEPT within the region (distant B)
    ##
    ##      3. Calculate interference B within the 0rec region
    ##
    ##      4. For regions with r between 0.1 and >0, compare result to normal calculateB and keep biggest result
    ##
    ##      5. Update B_values array


    ## linkage_block_max_rec_length = r * 100 # r * l
    ## linkage_block_min_selected_length = u * 10000


    ### THEN RUN
    interfering_L = 10000
    prior_B = 0.9

    b_values = calculateB_hri(interfering_L, prior_B)
    return
# ## To delete

    b_values = calculateB_linear(distance, gene_size)
    print(f"B for site {distance}bp away from {gene_size}bp region: {b_values}")

    if args.pop_change:
        if not quiet: print("Demographic change prior to B-calculation", b_values)
        b_values = get_Bcur(b_values)
        if not quiet: print("Demographic change applied to B-calculation", b_values)
    
    return