from Bvalcalc.core.calculateB import calculateB_linear, calculateB_recmap
from Bvalcalc.core.helpers.calc_R_len_dist import calc_R_lendist_for_chunks
import numpy as np

def calc_B_from_chunks(chunk_idx, chunk_size, chr_start, chr_size, num_chunks, 
                    precise_chunks, lperchunk, rec_rate_per_chunk, gc_rate_per_chunk, excluded_chunks = None):

    if excluded_chunks is not None:
        lperchunk_filtered = lperchunk.copy() # Branch lperchunk_filtered from lperchunk
        lperchunk_filtered[np.asarray(excluded_chunks, dtype=int)] = 0 # Set lperchunk to 0 for excluded chunks
    else: lperchunk_filtered = lperchunk # Else it's just lperchunk

    chunk_starts = chr_start + np.arange(num_chunks) * chunk_size
    chunk_ends = np.minimum(chunk_starts + chunk_size - 1, chr_size)
    chunk_mids = (chunk_ends + chunk_starts) / 2

    chunk_pseudoblockstart = chunk_mids - 0.5 * lperchunk_filtered # Pseudoblocks are the combined selected regions in each chunk, which are combined into a single central "pseudoblock" aka single selected region
    chunk_pseudoblockend = chunk_mids + 0.5 * lperchunk_filtered

    upstream_pseudochunk_mask = np.ones(num_chunks, dtype=bool)
    downstream_pseudochunk_mask = np.ones(num_chunks, dtype=bool)

    upstream_pseudochunk_mask[max(0, chunk_idx - precise_chunks):] = False # Mask for downstream blocks
    upstream_pseudochunk_mask[lperchunk_filtered==0] = False
    downstream_pseudochunk_mask[0:min(num_chunks, chunk_idx + precise_chunks + 1)] = False # Mask for upstream blocks
    downstream_pseudochunk_mask[lperchunk_filtered==0] = False
    
    relevant_upstream_psdc_lengths = lperchunk_filtered[upstream_pseudochunk_mask] # Use mask on pseudochunk lengths
    relevant_downstream_psdc_lengths = lperchunk_filtered[downstream_pseudochunk_mask] # Use mask on pseudochunk lengths

    relevant_upstream_pseudoblockends = chunk_pseudoblockend[upstream_pseudochunk_mask]
    relevant_downstream_pseudoblockstarts = chunk_pseudoblockstart[downstream_pseudochunk_mask]
    relevant_upstream_psdc_distances = chunk_mids[chunk_idx] - relevant_upstream_pseudoblockends - 1
    relevant_downstream_psdc_distances = relevant_downstream_pseudoblockstarts - chunk_mids[chunk_idx] - 1

    if rec_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE 
        # Get the indices for upstream and downstream pseudochunks
        chunk_rec_distances = (chunk_ends - chunk_starts + 1) * rec_rate_per_chunk
        upstream_indices = np.nonzero(upstream_pseudochunk_mask)[0]
        downstream_indices = np.nonzero(downstream_pseudochunk_mask)[0]

        upstream_rec_lengths, downstream_rec_lengths, upstream_rec_distances, downstream_rec_distances = calc_R_lendist_for_chunks(
            upstream_indices, downstream_indices, rec_rate_per_chunk, 
            relevant_upstream_psdc_lengths, relevant_downstream_psdc_lengths, 
            chunk_idx, chunk_size, relevant_upstream_pseudoblockends, relevant_downstream_pseudoblockstarts, 
            chunk_starts, chunk_ends, chunk_rec_distances, num_chunks
            ) # Get local r * lengths for length of, and distances to pseudoblocks for each chunk
        
    if gc_rate_per_chunk is not None: # IF GC_RATE MAP IS AVAILABLE 
        # Get the indices for upstream and downstream pseudochunks
        chunk_gc_distances = (chunk_ends - chunk_starts + 1) * gc_rate_per_chunk
        upstream_indices = np.nonzero(upstream_pseudochunk_mask)[0]
        downstream_indices = np.nonzero(downstream_pseudochunk_mask)[0]

        upstream_gc_lengths, downstream_gc_lengths, upstream_gc_distances, downstream_gc_distances = calc_R_lendist_for_chunks(
            upstream_indices, downstream_indices, gc_rate_per_chunk, 
            relevant_upstream_psdc_lengths, relevant_downstream_psdc_lengths, 
            chunk_idx, chunk_size, relevant_upstream_pseudoblockends, relevant_downstream_pseudoblockstarts, 
            chunk_starts, chunk_ends, chunk_gc_distances, num_chunks
            ) # Get local r * lengths for length of, and distances to pseudoblocks for each chunk


    #Run calculateBs
    if rec_rate_per_chunk is not None and gc_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE and GC IS AVAILABLE
        relevant_upstream_psdc_B = np.prod(calculateB_recmap(distance_to_element=relevant_upstream_psdc_distances, length_of_element=relevant_upstream_psdc_lengths, 
                                                             rec_distances=upstream_rec_distances, rec_lengths=upstream_rec_lengths, 
                                                               gc_distances=upstream_gc_distances, gc_lengths=upstream_gc_lengths))
        relevant_downstream_psdc_B = np.prod(calculateB_recmap(distance_to_element=relevant_downstream_psdc_distances, length_of_element=relevant_downstream_psdc_lengths, 
                                                               rec_distances=downstream_rec_distances, rec_lengths=downstream_rec_lengths, 
                                                               gc_distances=downstream_gc_distances, gc_lengths=downstream_gc_lengths))
    elif rec_rate_per_chunk is not None and gc_rate_per_chunk is None: # IF REC_RATE MAP IS AVAILABLE and GC NOT AVAILABLE
        relevant_upstream_psdc_B = np.prod(calculateB_recmap(distance_to_element=relevant_upstream_psdc_distances, length_of_element=relevant_upstream_psdc_lengths, 
                                                             rec_distances=upstream_rec_distances, rec_lengths=upstream_rec_lengths, 
                                                               gc_distances=None, gc_lengths=None))
        relevant_downstream_psdc_B = np.prod(calculateB_recmap(distance_to_element=relevant_downstream_psdc_distances, length_of_element=relevant_downstream_psdc_lengths, 
                                                               rec_distances=downstream_rec_distances, rec_lengths=downstream_rec_lengths, 
                                                               gc_distances=None, gc_lengths=None))
    elif rec_rate_per_chunk is None and gc_rate_per_chunk is not None: # IF REC_RATE MAP NOT AVAILABLE and GC IS AVAILALBE
        relevant_upstream_psdc_B = np.prod(calculateB_recmap(distance_to_element=relevant_upstream_psdc_distances, length_of_element=relevant_upstream_psdc_lengths, 
                                                             rec_distances=None, rec_lengths=None, 
                                                               gc_distances=upstream_gc_distances, gc_lengths=upstream_gc_lengths))
        relevant_downstream_psdc_B = np.prod(calculateB_recmap(distance_to_element=relevant_downstream_psdc_distances, length_of_element=relevant_downstream_psdc_lengths, 
                                                               rec_distances=None, rec_lengths=None, 
                                                               gc_distances=downstream_gc_distances, gc_lengths=downstream_gc_lengths))
    else: # NEITHER REC_MAP NOR GC_MAP AVAILABLE
        relevant_upstream_psdc_B = np.prod(calculateB_linear(relevant_upstream_psdc_distances, relevant_upstream_psdc_lengths))
        relevant_downstream_psdc_B = np.prod(calculateB_linear(relevant_downstream_psdc_distances, relevant_downstream_psdc_lengths))

    return relevant_downstream_psdc_B * relevant_upstream_psdc_B # Return a single B value that applies to all sites in focal chunk