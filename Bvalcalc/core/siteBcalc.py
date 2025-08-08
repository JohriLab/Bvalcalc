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

    ##      0. Make a normal B-map.

    ##      1. Find consecutive chunks with local rec_rate below 0.1 * r with int_L > 0 to have interference calculation applied within them.
    import numpy as np
    r = 1e-8
    r_hri_threshold = 0.1 # only run HRI B when chunk r is less than: r * r_hri_threshold

    chunk_rec_rate = np.array([0.001 * r, 0.5 * r])
    print(chunk_rec_rate[chunk_rec_rate > r_hri_threshold * r], 'Hi')

    import sys
    sys.exit()

    ##      1b. Combine the single chunks with adjacent chunks   
    ## 
    ##      2. Re-calculate B from everything EXCEPT within the region (distant B)
    ##
    ##      3. Calculate interference B within the 0rec region
    ##
    ##      4. For regions with r between 0.1 and >0, compare result to normal calculateB and keep biggest result
    ##
    ##      5. Update B_values array


    ##  For sites outside the interference region, calculating B from sites within the interference region
    ##
    ##  Option 1. Use the B from within the interference region and modify it by r accordingly
    ##
    ##  Option 2. Ignore it
    ## 
    ##  Option 3a. Calculate interference B for all interfering regions separately.
    ##  Option 3b. Calculate B as normal for all non-interfering regions
    ##  Option 3c. For within interference region, use interference B
    ##  Option 3d. For outside interference region, use interference B modified by rec length and rec distance
    ##
    ##  Option 4. Calculate b from interfering region as if there's no interference, but algorithmically ensure there's no local dip at the edge of non-recombining regions.
    ##
    ##  Consequences: BGS from interference region is overestimated, as it's not accounting for HRI reducing efficacy of selection in interference region
    ##                  Response: For this reason, our results will not be accurate if a large proportion of the genome experiencing selection also experiences very low recombination.
    ##                  Response: In other cases, this is only expected to substantially affect nearby flanking regions to a HRI region, particularly due to the effects of weakly and moderately deleterious mutations (f1,f2 class).
    ##                  Response: To illustrate this effect, consider a HRI region which has a recombination rate of 0 and 10kb of selected sites. 
    ##                             B' can be well approximated for these HRI region sites by Bvalcalc, even taking into account an accurate estimate of N0 from classic BGS effects (B) across the genome. 
    ##                             However, within the HRI region with r = 0, the recombinant distance is 0. 
    ##                             Which means for the first site outside the HRI region, B will be calculated as though there is 10kb of selected sites, but all of them are 1bp away.
    ##                             This is obviously ineffective, a more rigorous approach would calculate B for these sites considering the lessened effects of selection in the region due to HRI. 
    ##                             On a deeper level, this is an issue that arises due to the shortcut of partitioning interference regions into discrete categories in spite of the reality that the borders of regions experiencing HRI are diffuse.
    ##                             Nevertheless, while Bvalcalc does not implement an exact solution to this quandry, we have implemented a practical algorithm to smooth the borders between HRI regions and immediately adjacent sites.
    ##                             This algorithm runs after B is calculating for all non-HRI sites considering classic BGS from all sites, and after B' is calculated for each HRI region. 
    ##                                  For each site j distance from the HRI region, where j starts at 1 (bp), compare the B of site j to B' of the flanking HRI region. 
    #                                   If B' is higher than B, replace B with B' for this site and increment j. 
    #                                   Continue replacing B with B' until you reach a site distance j away for which B is higher than B', which should happen rapidly once the local effect from f1 and f2 mutations is diminished by recombination
    #                                   To avoid edge cases where B decreases due to proximity to a new selected region placed a few bp distance from the HRI region, break also if B at site j is higher than at site j+1 (as in, B is decreasing with increasing distance from the HRI region)
    ##
    ##  4a. Calculate 
    ##
    ##


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