import numpy as np


def calc_R_lengths(blockstart, blockend, rec_rate_per_chunk, calc_start, calc_end, chunk_size, chunk_num):
    """
    Calculates the weighted lengths of each conserved block (gene), so that for example if the mean 
    recombination rate across the block is 0.5, this will return the length of the block multiplied by 0.5
    """
    num_chunks = (calc_end - calc_start) // chunk_size
    # Build chunk boundaries (note: length = num_chunks + 1)
    chunk_starts = calc_start + np.arange(0, num_chunks + 1) * chunk_size
    chunk_left  = chunk_starts           # shape: (num_chunks+1,)
    chunk_right = chunk_starts + chunk_size  # shape: (num_chunks+1,)

    # Ensure block boundaries are numpy arrays (and 1D)
    blockstart = np.asarray(blockstart).flatten()
    blockend   = np.asarray(blockend).flatten()

    # Compute the overlap between each block and each chunk interval:
    # For block i and chunk j, the overlap is:
    #   max(0, min(blockend[i], chunk_right[j]) - max(blockstart[i], chunk_left[j]))
    overlap = np.maximum(0, np.minimum(blockend[:, None], chunk_right[None, :]) -
                           np.maximum(blockstart[:, None], chunk_left[None, :]))
    weighted_overlap = overlap * rec_rate_per_chunk[None, :] # Multiply by the recombination rate for each chunk
    weighted_sum = np.sum(weighted_overlap, axis=1) # Sum over the chunk intervals for each block
    
    return weighted_sum

def calc_R_distances(precise_blockstart, precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, pos_chunk_clean, chunk_num, chunk_start):

    num_chunks = (precise_region_end - precise_region_start) // chunk_size
    chunk_starts = precise_region_start + np.arange(0, num_chunks + 1) * chunk_size
    chunk_ends = np.minimum(chunk_starts + chunk_size, precise_region_end)
    this_chunk_idx = np.where(chunk_starts == chunk_start)[0][0] # The ID of this chunk in the chunk_starts array, e.g. if precise_chunks = 3, this will be [3] for chunk_num > 2
    chunk_end = chunk_ends[this_chunk_idx] # Focal chunk's end
    blockstart_chunks = (precise_blockstart - precise_region_start) // chunk_size
    blockend_chunks = (precise_blockend - precise_region_start) // chunk_size
    blockend_rec_distances = np.empty((len(blockend_chunks), len(pos_chunk_clean)))
    blockstart_rec_distances = np.empty((len(blockstart_chunks), len(pos_chunk_clean)))

## FOR LOOP COULD BE IMPORVED BY USING NP ARRAY OPERATIONS INSTEAD!!!

    for block_idx in range(len(blockend_chunks)): # For blockends (i.e. upstream genes)
        rec_distance_overlapped, rec_distance_blockchunk = 0, 0 # Set to 0 to allow for sums even when not relevant
        inchunk_distances = np.minimum(pos_chunk_clean - precise_blockend[block_idx], 
                                       pos_chunk_clean - chunk_start) # To blockend if block is within same chunk, else to chunk start
        rec_distance_focalchunk = inchunk_distances * precise_rates[this_chunk_idx]
        isin_diffchunk = (blockend_chunks[block_idx] < this_chunk_idx) # 0 if in same chunk, 1 if in different chunk

        chunk_edge_distances = chunk_ends[blockend_chunks[block_idx]] - precise_blockend[block_idx]  # Distances from the end of the block to the end of its chunk
        rec_distance_blockchunk = chunk_edge_distances * precise_rates[blockend_chunks[block_idx]] # Rec_distance to end of block's chunk NOT spanned chunks

        overlapped_chunks = np.arange(blockend_chunks[block_idx] + 1, this_chunk_idx)
        rec_distance_overlapped = np.sum(precise_rates[overlapped_chunks] * chunk_size) # Rec_distance in chunks that are overlapped

        total_rec_distances = np.array(rec_distance_focalchunk + isin_diffchunk * (rec_distance_blockchunk + rec_distance_overlapped))
                    
        blockend_rec_distances[block_idx] = total_rec_distances

    for block_idx in range(len(blockstart_chunks)): # For blockstarts (i.e. downstream genes)
        rec_distance_overlapped, rec_distance_blockchunk = 0, 0 # Set to 0 to allow for sums even when not relevant
        inchunk_distances = np.minimum(precise_blockstart[block_idx] - pos_chunk_clean, 
                                       chunk_end - pos_chunk_clean) # To blockstart if block is within same chunk, else to chunk start
        rec_distance_focalchunk = inchunk_distances * precise_rates[this_chunk_idx]
        isin_diffchunk = (blockstart_chunks[block_idx] > this_chunk_idx) # 0 if in same chunk, 1 if in different chunk

        chunk_edge_distances = precise_blockstart[block_idx] - chunk_starts[blockstart_chunks[block_idx]]   # Distances from the start of the block to the start of its chunk
        rec_distance_blockchunk = chunk_edge_distances * precise_rates[blockstart_chunks[block_idx]] # Rec_distance to start of block's chunk NOT spanned chunks

        overlapped_chunks = np.arange(this_chunk_idx + 1, blockstart_chunks[block_idx])
        rec_distance_overlapped = np.sum(precise_rates[overlapped_chunks] * chunk_size) # Rec_distance in chunks that are overlapped

        total_rec_distances = np.array(rec_distance_focalchunk + isin_diffchunk * (rec_distance_blockchunk + rec_distance_overlapped))

        blockstart_rec_distances[block_idx] = total_rec_distances

    return blockend_rec_distances, blockstart_rec_distances

def calc_R_lendist_for_chunks(upstream_indices, downstream_indices, rec_rate_per_chunk, relevant_upstream_psdc_lengths, relevant_downstream_psdc_lengths, chunk_index, chunk_size, relevant_upstream_pseudoblockends, relevant_downstream_pseudoblockstarts, chunk_starts, chunk_ends, chunk_rec_distances, num_chunks):

    ## Calculate relevant upstream and downstream rec lengths of pseudoblocks
    upstream_rec_rates = rec_rate_per_chunk[upstream_indices] # Relevant rec rates for pseudochunks upstream
    upstream_rec_lengths = upstream_rec_rates * relevant_upstream_psdc_lengths
    downstream_rec_rates = rec_rate_per_chunk[downstream_indices] # Relevant rec rates for pseudochunks downstream
    downstream_rec_lengths = downstream_rec_rates * relevant_downstream_psdc_lengths

    ## Calculate relevant upstream rec distances!
    ## CURRENTLY TAKING A LONG TIME TO CALCULATE, could potentially improve with array functions
    mean_rec_distance_focalchunk = rec_rate_per_chunk[chunk_index] * chunk_size / 2 - 0.5 # Note that this is distance to middle of focal chunk.

    if chunk_index == num_chunks -1: # If it's the final chunk (which may be not be full chunk_size length)
        end_focalchunk_distance = (chunk_ends[chunk_index] - chunk_starts[chunk_index])/2 # I think needs to be shifted by 1bp
        mean_rec_distance_focalchunk = end_focalchunk_distance

    upstream_distance_blockchunk = chunk_ends[upstream_indices] - relevant_upstream_pseudoblockends
    upstream_rec_distance_blockchunk = upstream_distance_blockchunk * rec_rate_per_chunk[upstream_indices] # This is rec distance from edge of pseudoblock to its chunk end

    upstream_overlapped_indices = [np.arange(u + 1, chunk_index) for u in upstream_indices]
    upstream_overlapped_rec_distances = np.array([chunk_rec_distances[idx].sum() for idx in upstream_overlapped_indices]) # This is rec distance spanned in fully overlapped chunks

    upstream_rec_distances = mean_rec_distance_focalchunk + upstream_rec_distance_blockchunk + upstream_overlapped_rec_distances # Combined rec distance from middle of focal chunk to edge of pseudo"blocks" upstream

    ## Calculate downstream rec distances!
    ## CURRENTLY TAKING A LONG TIME TO CALCULATE, could potentially improve with array functions
    downstream_distance_blockchunk = relevant_downstream_pseudoblockstarts - chunk_starts[downstream_indices]
    downstream_rec_distance_blockchunk = downstream_distance_blockchunk * rec_rate_per_chunk[downstream_indices] # This is rec distance from edge of pseudoblock to its chunk start
    
    downstream_overlapped_indices = [np.arange(chunk_index + 1, d) for d in downstream_indices]
    downstream_overlapped_rec_distances = np.array([chunk_rec_distances[idx].sum() for idx in downstream_overlapped_indices]) # This is rec distance spanned in fully overlapped chunks

    downstream_rec_distances = mean_rec_distance_focalchunk + downstream_rec_distance_blockchunk + downstream_overlapped_rec_distances # Combined rec distance from middle of focal chunk to edge of pseudo"blocks" upstream


    return upstream_rec_lengths, downstream_rec_lengths, upstream_rec_distances, downstream_rec_distances


## Can remove??
# def calc_R_lendist_for_genes(this_chunk_blockstart, this_chunk_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, chunk_num, chunk_start, precise_chunks):
#     num_chunks = (precise_region_end - precise_region_start) // chunk_size
#     chunk_starts = precise_region_start + np.arange(0, num_chunks + 1) * chunk_size
#     chunk_ends = np.minimum(chunk_starts + chunk_size, precise_region_end)
#     this_chunk_idx = np.where(chunk_starts == chunk_start)[0][0] # The ID of this chunk in the chunk_starts array, e.g. if precise_chunks = 3, this will be [3] for chunk_num > 2
#     right_rec_distances = 1 * precise_rates[this_chunk_idx]
#     left_rec_distances = 1 * precise_rates[this_chunk_idx]

    
#     print("in recmaphandler", precise_rates)
#     right_rec_lengths = 3

#     return right_rec_lengths, right_rec_distances, left_rec_distances