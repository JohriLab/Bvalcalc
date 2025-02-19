from helperScripts.calculateB import calculateB
import numpy as np

def calcBFromChunks(chunk_index, chunk_size, blockstart, blockend, chr_start, chr_end, num_chunks, precise_chunks, lperchunk):

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

    relevant_upstream_psdc_B = np.prod(calculateB(relevant_upstream_psdc_distances, relevant_upstream_psdc_lengths, rdistance_to_element=None, rlength_of_element=None))
    relevant_downstream_psdc_B = calculateB(relevant_downstream_psdc_distances, relevant_downstream_psdc_lengths, rdistance_to_element=None, rlength_of_element=None)

    return np.prod(relevant_downstream_psdc_B) * relevant_upstream_psdc_B # Return B that applies to all sites in a chunk