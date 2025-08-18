import numpy as np

def calc_R_lengths(blockstart, blockend, rec_rate_per_chunk, calc_start, calc_end, chunk_size):
    """
    Calculates the weighted lengths of each conserved block (gene), so that for example if the mean 
    recombination rate across the block is 0.5, this will return the length of the block multiplied by 0.5
    """
    num_precise_chunks = (calc_end - calc_start - 1) // chunk_size # -1 to fix extra chunk bug. May need to revert.

    # Build chunk boundaries (note: length = num_chunks + 1)
    chunk_starts = calc_start + np.arange(0, num_precise_chunks + 1) * chunk_size
    next_chunk_starts = chunk_starts + chunk_size  # shape: (num_chunks+1,)

    # Compute the overlap between each block and each chunk interval:
    # For block i and chunk j, the overlap is:
    #   max(0, min(blockend[i], chunk_right[j]) - max(blockstart[i], chunk_left[j]))
    overlap = np.maximum(0, np.minimum(blockend[:, None], next_chunk_starts[None, :]) -
                           np.maximum(blockstart[:, None], chunk_starts[None, :]))
    weighted_overlap = overlap * rec_rate_per_chunk[None, :] # Multiply by the recombination rate for each chunk
    weighted_sum = np.sum(weighted_overlap, axis=1) # Sum over the chunk intervals for each block
    
    return weighted_sum

def calc_R_distances(
    precise_blockstart, precise_blockend, precise_rates,
    precise_region_start, precise_region_end, chunk_size,
    pos_chunk_clean, chunk_start
):
    

    num_precise_chunks = (precise_region_end - precise_region_start) // chunk_size
    chunk_starts = precise_region_start + np.arange(0, num_precise_chunks + 1) * chunk_size
    chunk_ends = np.minimum(chunk_starts + chunk_size, precise_region_end)
    this_chunk_idx = np.where(chunk_starts == chunk_start)[0][0]
    chunk_end = chunk_ends[this_chunk_idx]

    blockstart_chunks = (precise_blockstart - precise_region_start) // chunk_size
    blockend_chunks = (precise_blockend - precise_region_start) // chunk_size

    pos_broadcast = pos_chunk_clean[None, :]  # shape (1, P)
    # BLOCKEND (upstream) distances
    inchunk_distances_end = np.minimum(
        pos_broadcast - precise_blockend[:, None],
        pos_broadcast - chunk_start
    )
    rec_distance_focal_end = inchunk_distances_end * precise_rates[this_chunk_idx]

    # Mask out focal distances if block is in a future chunk
    same_or_earlier = (blockend_chunks <= this_chunk_idx)[:, None]
    rec_distance_focal_end *= same_or_earlier

    is_diffchunk_end = (blockend_chunks < this_chunk_idx).astype(float)[:, None]
    chunk_edge_dist_end = chunk_ends[blockend_chunks] - precise_blockend
    rec_blockchunk_end = chunk_edge_dist_end * precise_rates[blockend_chunks]
    blockend_overlap_dists = np.array([
        np.sum(precise_rates[blockend_chunks[i]+1:this_chunk_idx] * chunk_size)
        for i in range(len(blockend_chunks))
    ])
    blockend_rec_distances = rec_distance_focal_end + is_diffchunk_end * (
        rec_blockchunk_end[:, None] + blockend_overlap_dists[:, None]
    )

    # BLOCKSTART (downstream) distances
    inchunk_distances_start = np.minimum(
        precise_blockstart[:, None] - pos_broadcast,
        chunk_end - pos_broadcast
    )
    rec_distance_focal_start = inchunk_distances_start * precise_rates[this_chunk_idx]

    # Mask out focal distances if block is in a previous chunk
    same_or_later = (blockstart_chunks >= this_chunk_idx)[:, None]
    rec_distance_focal_start *= same_or_later

    is_diffchunk_start = (blockstart_chunks > this_chunk_idx).astype(float)[:, None]
    chunk_edge_dist_start = precise_blockstart - chunk_starts[blockstart_chunks]
    rec_blockchunk_start = chunk_edge_dist_start * precise_rates[blockstart_chunks]
    blockstart_overlap_dists = np.array([
        np.sum(precise_rates[this_chunk_idx+1:blockstart_chunks[i]] * chunk_size)
        for i in range(len(blockstart_chunks))
    ])
    blockstart_rec_distances = rec_distance_focal_start + is_diffchunk_start * (
        rec_blockchunk_start[:, None] + blockstart_overlap_dists[:, None]
    )

    return blockend_rec_distances, blockstart_rec_distances

def calc_R_lendist_for_chunks(upstream_indices, downstream_indices, rec_rate_per_chunk, relevant_upstream_psdc_lengths, relevant_downstream_psdc_lengths, chunk_idx, chunk_size, relevant_upstream_pseudoblockends, relevant_downstream_pseudoblockstarts, chunk_starts, chunk_ends, chunk_rec_distances, num_chunks):

    ## Calculate relevant upstream and downstream rec lengths of pseudoblocks
    upstream_rec_rates = rec_rate_per_chunk[upstream_indices] # Relevant rec rates for pseudochunks upstream
    upstream_rec_lengths = upstream_rec_rates * relevant_upstream_psdc_lengths
    downstream_rec_rates = rec_rate_per_chunk[downstream_indices] # Relevant rec rates for pseudochunks downstream
    downstream_rec_lengths = downstream_rec_rates * relevant_downstream_psdc_lengths

    ## Calculate relevant upstream rec distances!
    mean_rec_distance_focalchunk = rec_rate_per_chunk[chunk_idx] * chunk_size / 2 - 0.5 # Note that this is distance to middle of focal chunk.

    if chunk_idx == num_chunks -1: # If it's the final chunk (which may be not be full chunk_size length)
        end_focalchunk_distance = (chunk_ends[chunk_idx] - chunk_starts[chunk_idx])/2 # I think needs to be shifted by 1bp
        mean_rec_distance_focalchunk = end_focalchunk_distance

    upstream_distance_blockchunk = chunk_ends[upstream_indices] - relevant_upstream_pseudoblockends
    upstream_rec_distance_blockchunk = upstream_distance_blockchunk * rec_rate_per_chunk[upstream_indices] # This is rec distance from edge of pseudoblock to its chunk end

    chunk_cumsum = np.concatenate([[0], np.cumsum(chunk_rec_distances)])
    upstream_start = np.array(upstream_indices) + 1
    upstream_end = chunk_idx
    upstream_overlapped_rec_distances = chunk_cumsum[upstream_end] - chunk_cumsum[upstream_start]  # This is rec distance spanned in fully overlapped chunks

    upstream_rec_distances = mean_rec_distance_focalchunk + upstream_rec_distance_blockchunk + upstream_overlapped_rec_distances # Combined rec distance from middle of focal chunk to edge of pseudo"blocks" upstream

    ## Calculate downstream rec distances!
    downstream_distance_blockchunk = relevant_downstream_pseudoblockstarts - chunk_starts[downstream_indices]
    downstream_rec_distance_blockchunk = downstream_distance_blockchunk * rec_rate_per_chunk[downstream_indices] # This is rec distance from edge of pseudoblock to its chunk start
    
    downstream_start = chunk_idx + 1
    downstream_end = np.array(downstream_indices)
    downstream_overlapped_rec_distances = chunk_cumsum[downstream_end] - chunk_cumsum[downstream_start] # This is rec distance spanned in fully overlapped chunks

    downstream_rec_distances = mean_rec_distance_focalchunk + downstream_rec_distance_blockchunk + downstream_overlapped_rec_distances # Combined rec distance from middle of focal chunk to edge of pseudo"blocks" upstream

    return upstream_rec_lengths, downstream_rec_lengths, upstream_rec_distances, downstream_rec_distances