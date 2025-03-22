import numpy as np

def calcBInGenes(chunk_index, chunk_size, chr_start, chr_end, num_chunks, 
                    precise_chunks, lperchunk, gene_mask, rec_rate_per_chunk, gc_rate_per_chunk):
    if chunk_index == 3:
        print("gottem", np.sum(gene_mask, axis=1))

        # NOW NEED TO CALCULATE B FOR WITHIN THE GENE!!