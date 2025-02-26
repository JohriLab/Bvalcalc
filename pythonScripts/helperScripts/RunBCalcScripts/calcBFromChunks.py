from helperScripts.calculateB import calculateB_linear, calculateB_recmap
import numpy as np

def calcBFromChunks(chunk_index, chunk_size, blockstart, blockend, chr_start, chr_end, num_chunks, precise_chunks, lperchunk, rec_rate_per_chunk):


    chunk_mids = chr_start + (np.arange(num_chunks) + 0.5) * chunk_size
    chunk_pseudoblockstart = chunk_mids - 0.5 * lperchunk
    chunk_pseudoblockend = chunk_mids + 0.5 * lperchunk

    upstream_pseudochunk_mask = np.ones(num_chunks, dtype=bool)
    downstream_pseudochunk_mask = np.ones(num_chunks, dtype=bool)

    upstream_pseudochunk_mask[max(0, chunk_index - precise_chunks):] = False # Mask for downstream blocks
    upstream_pseudochunk_mask[lperchunk==0] = False
    downstream_pseudochunk_mask[0:min(num_chunks, chunk_index + precise_chunks + 1)] = False # Mask for upstream blocks
    downstream_pseudochunk_mask[lperchunk==0] = False
    
    relevant_upstream_psdc_lengths = lperchunk[upstream_pseudochunk_mask] # Use mask on pseudochunk lengths
    relevant_downstream_psdc_lengths = lperchunk[downstream_pseudochunk_mask] # Use mask on pseudochunk lengths

    relevant_upstream_pseudoblockends = chunk_pseudoblockend[upstream_pseudochunk_mask]
    relevant_downstream_pseudoblockstarts = chunk_pseudoblockstart[downstream_pseudochunk_mask]
    relevant_upstream_psdc_distances = chunk_mids[chunk_index] - relevant_upstream_pseudoblockends
    relevant_downstream_psdc_distances = relevant_downstream_pseudoblockstarts - chunk_mids[chunk_index]

    if rec_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE 
        # Get the indices for upstream and downstream pseudochunks
        upstream_indices = np.nonzero(upstream_pseudochunk_mask)[0]
        downstream_indices = np.nonzero(downstream_pseudochunk_mask)[0]

        upstream_rec_rates = rec_rate_per_chunk[upstream_indices] # Relevant rec rates for pseudochunks upstream
        downstream_rec_rates = rec_rate_per_chunk[downstream_indices] # Relevant rec rates for pseudochunks downstream

        relevant_upstream_psdc_B = np.prod(calculateB_recmap(relevant_upstream_psdc_distances, relevant_upstream_psdc_lengths, rec_gene_modifier=upstream_rec_rates))
        relevant_downstream_psdc_B = calculateB_recmap(relevant_downstream_psdc_distances, relevant_downstream_psdc_lengths, rec_gene_modifier=downstream_rec_rates)
    else:
        relevant_upstream_psdc_B = np.prod(calculateB_linear(relevant_upstream_psdc_distances, relevant_upstream_psdc_lengths))
        relevant_downstream_psdc_B = calculateB_linear(relevant_downstream_psdc_distances, relevant_downstream_psdc_lengths)
        


    return np.prod(relevant_downstream_psdc_B) * relevant_upstream_psdc_B # Return B that applies to all sites in a chunk