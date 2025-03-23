import numpy as np
from helperScripts.calculateB import calculateB_linear

def calcBInGenes(chunk_index, chunk_size, chr_start, chr_end, num_chunks, 
                    precise_chunks, lperchunk, gene_mask, rec_rate_per_chunk, gc_rate_per_chunk):
    



    # For each gene in this chunk:
    #
    # # Make a left gene array
    # # # Run Bcalc on the left gene
    # # # # Save in an array with position for each B
    #
    # # Make a right gene array
    # # # Run Bcalc on the right gene
    # # # # Save in an array with position for each B
    #

    length_of_genes_in_chunk = 9000
    # np.arange(1,000)

    if chunk_index == 3:
        # print("hi6", np.sum(gene_mask))
        # print("left_gene_distances",np.arange(1,1000))
        # calculateB_linear(distance, length )
        print("gottem", np.sum(gene_mask[4], axis=1))

        # NOW NEED TO CALCULATE B FOR WITHIN THE GENE!!