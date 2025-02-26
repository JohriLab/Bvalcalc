import sys
import numpy as np
from helperScripts.calculateB import calculateB_linear, calculateB_recmap
from helperScripts.RunBCalcScripts.calcBFromChunks import calcBFromChunks
from helperScripts.RunBCalcScripts.recmapHandler import calcRLengths
from helperScripts.RunBCalcScripts.recmapHandler import calcRDistances

def process_single_chunk(chunk_num, chunk_size, blockstart, blockend,
                         chr_start, chr_end, num_chunks, precise_chunks,
                         lperchunk, b_values, rec_rate_per_chunk=None):
    
    # Compute chunk boundaries
    chunk_start = chr_start + chunk_num * chunk_size
    chunk_end   = min(chunk_start + chunk_size, chr_end)
    
    # Array of positions in this chunk
    pos_chunk = np.arange(chunk_start, chunk_end)

    # Slice b_values for this chunk
    start_idx = chunk_start - chr_start
    end_idx   = chunk_end - chr_start
    chunk_slice = b_values[start_idx:end_idx]
    # chunk_slice.shape == (len(pos_chunk),)

    # Make a mask for positions that are NOT NaN
    not_nan_mask = ~np.isnan(chunk_slice)
    
    if not np.any(not_nan_mask):
        # Everything is already NaN in this chunk, no need to do anything
        print(f"No neutral sites in chunk {chunk_num}: {chunk_start}-{chunk_end}")
        return b_values
    
    # Filter positions to skip pre-labeled NaNs
    pos_chunk_clean   = pos_chunk[not_nan_mask] # Positions of neutral sites
    chunk_slice_clean = chunk_slice[not_nan_mask]

    # == 1) If you rely on B_from_distant_chunks, compute that as usual ==
    B_from_distant_chunks = calcBFromChunks(
        chunk_num, chunk_size,
        blockstart, blockend,
        chr_start, chr_end,
        num_chunks, precise_chunks,
        lperchunk,
        rec_rate_per_chunk
    )

    # == 2) Identify blocks in the "precise region" (unchanged) ==
    precise_region_start = np.maximum(chr_start, chr_start + (chunk_num - precise_chunks) * chunk_size)
    precise_region_end   = np.minimum(chr_end, chr_start + (chunk_num + 1 + precise_chunks) * chunk_size - 1)

    precise_blockregion_mask = (
        (precise_region_end   >= blockstart) &
        (precise_region_start <= blockend)
    )
    precise_blockstart = np.clip(blockstart[precise_blockregion_mask],
                                 a_min=precise_region_start, a_max=precise_region_end)
    precise_blockend   = np.clip(blockend[precise_blockregion_mask],
                                 a_min=precise_region_start, a_max=precise_region_end)
    # print(precise_blockend)



    # # == 3) Do distance calculations ONLY for non-NaN sites ==
    downstream_mask = (pos_chunk_clean < precise_blockstart[:, None])
    upstream_mask   = (pos_chunk_clean > precise_blockend[:, None])
    flanking_mask   = downstream_mask | upstream_mask

    physical_distances_upstream   = pos_chunk_clean[None, :] - precise_blockend[:, None]
    if chunk_num == 12:
        print("HERE1st", np.sum(physical_distances_upstream), chunk_num)
    physical_distances_downstream = precise_blockstart[:, None] - pos_chunk_clean[None, :]
    physical_distances = np.where(
        flanking_mask,
        np.where(upstream_mask, physical_distances_upstream, physical_distances_downstream),
        np.nan
    )

    flat_distances = physical_distances[flanking_mask]

    physical_lengths = precise_blockend - precise_blockstart
    flat_lengths   = np.repeat(physical_lengths, flanking_mask.sum(axis=1))
    nonzero_mask = flat_lengths != 0 # Remove genes of length 0
    flat_distances = flat_distances[nonzero_mask]
    flat_lengths   = flat_lengths[nonzero_mask]
    
    if rec_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE 
        precise_rates = rec_rate_per_chunk[np.maximum(0, chunk_num - precise_chunks):np.minimum(num_chunks, chunk_num + precise_chunks + 1)]
        rec_lengths = calcRLengths(precise_blockstart, precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, chunk_num)
        rec_distances_upstream, rec_distances_downstream = calcRDistances(precise_blockstart, precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, pos_chunk_clean, chunk_num, chunk_start)
        if chunk_num == 12:
            print("HERE", np.sum(rec_distances_upstream), chunk_num)
        rec_distances = np.where(
            flanking_mask,
            np.where(upstream_mask, rec_distances_upstream, rec_distances_downstream),
            np.nan
        )
        flat_rec_distances = rec_distances[flanking_mask]
        flat_rec_lengths   = np.repeat(rec_lengths, flanking_mask.sum(axis=1))
        nonzero_rec_mask = flat_rec_lengths != 0 # Remove genes of length 0
        flat_rec_distances = flat_rec_distances[nonzero_rec_mask]
        flat_rec_lengths   = flat_rec_lengths[nonzero_rec_mask]

    if rec_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE 
        flank_B = calculateB_recmap(flat_distances, flat_lengths, flat_rec_distances, flat_rec_lengths)
        # print(np.sum(flat_distances), np.sum(flat_lengths), np.sum(flat_rec_distances), np.sum(flat_rec_lengths))
    else:
        flank_B = calculateB_linear(flat_distances, flat_lengths)

    # If there are no elements left, default flank_B to 1
    if flat_distances.size == 0 or flat_lengths.size == 0:
        flank_B = 1
        print("No elements left 'process_single_chunk'")

    true_indices = np.where(flanking_mask)
    unique_indices, inverse_indices = np.unique(true_indices[1], return_inverse=True)

    aggregated_B = np.ones_like(unique_indices, dtype=np.float64)
    np.multiply.at(aggregated_B, inverse_indices, flank_B)

    # == 4) Write results only to the non-NaN slice ==
    global_indices_clean = pos_chunk_clean[unique_indices] - chr_start
    # Now map them into chunk_slice_clean indices:
    chunk_slice_clean_indices = unique_indices  # because pos_chunk_clean -> chunk_slice_clean is 1:1
    chunk_slice_clean[chunk_slice_clean_indices] *= (aggregated_B * B_from_distant_chunks)

    # Put the updated (non-NaN) slice back into the original b_values
    # The NaNs remain untouched
    chunk_slice[not_nan_mask] = chunk_slice_clean
    mean_chunk_b = np.nanmean(chunk_slice) # Mean B for chunk

    # print(f"Processing chunk: {pos_chunk.min()} - {pos_chunk.max()}")
    #         # Check if recombination rate data is provided and print it for the current chunk.
    # if rec_rate_per_chunk is not None:
    #     rec_rate = rec_rate_per_chunk[chunk_num]
    #     print(f"Chunk {chunk_num}: recombination rate = {rec_rate}")
    # print(f"B from distant chunks: {B_from_distant_chunks}")
    # print(f"Number of relevant genes: {len(precise_blockstart)}")
    # # print(f"Relevant blocks: {precise_blockstart}, {precise_blockend}")
    # print(f"Number of neutral sites in chunk [{chunk_start}-{chunk_end}): {np.isnan(chunk_slice).sum()}")
    # print(f"Aggregated B values for chunk: {aggregated_B}")
    # print(f"Mean B value for chunk {chunk_num}: [{chunk_start}-{chunk_end}]: {mean_chunk_b}")

    return b_values