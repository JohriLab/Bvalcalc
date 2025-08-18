from Bvalcalc.core.helpers.demography_helpers import get_Bcur
from Bvalcalc.utils.bin_outputs import bin_outputs
import numpy as np
import os

def geneBcalc(args):    
    element_size, flank_len, quiet = args.element_size, args.flank_len, args.quiet
    
    import Bvalcalc.utils.dfe_helper as dfe_helper
    dfe_helper.CONSTANT_DFE = args.constant_dfe # Update DFE if --constant_dfe
    dfe_helper.GAMMA_DFE = args.gamma_dfe # Update DFE if --gamma_dfe
    from Bvalcalc.core.calculateB import calculateB_linear # Import calculateB with updated gamma DFE if needed

    if not quiet: 
        print(f"====== P A R A M E T E R S =========================")
        print(f"Distribution of fitness effects (DFE): {flank_len}bp")
        print(f"Length of element under selection: {element_size}bp")
        print(f"Length of flanking neutral region: {flank_len}bp")

    print(f"====== S T A R T I N G ===== C A L C ===============")
    b_values = calculateB_linear(np.arange(1, flank_len, 1, dtype = int), element_size) # B for flank region

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
    
    ## !!The following (prior to return) calculates output_data WITH the sites in the gene from the edge to the middle of the gene
    ## !!Not currently supported because it'll require major updates to plotting --gene, worth adding when plotB is refactored

        # distance_to_middle_of_gene = int((element_size)/2)
        # left_gene_lengths = np.arange(distance_to_middle_of_gene, distance_to_middle_of_gene*2, 1, dtype = int)
        # right_gene_lengths = np.arange(1, distance_to_middle_of_gene + 1, 1, dtype = int)[::-1] #Flip array at the end so it starts with center point

        # right_gene_b_values = calculateB_linear(distance_to_element = 0, length_of_element = right_gene_lengths) 
        # ##Distance = 0 includes BGS from focal site, doesn't make logical sense but avoids weird increase in B as selected elements have 1 less site when calculating within them
        # left_gene_b_values = calculateB_linear(distance_to_element = 1, length_of_element = left_gene_lengths)
        
        # in_gene_b_values = right_gene_b_values*left_gene_b_values
        # print(element_size, distance_to_middle_of_gene, np.shape(left_gene_b_values), np.shape(right_gene_b_values), left_gene_lengths, right_gene_lengths, in_gene_b_values)

        # combined_b_values = np.append(in_gene_b_values, b_values)
        # print(combined_b_values, np.shape(combined_b_values))
        # with_gene_output_data = np.column_stack((np.arange(1 - distance_to_middle_of_gene, len(combined_b_values) + 1 - distance_to_middle_of_gene, 1, dtype = int), combined_b_values))

        # print("Kaizo", output_data, with_gene_output_data)

    return output_data
