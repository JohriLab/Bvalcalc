from helperScripts.calculateB import calculateB_linear, calculateB_recmap
from helperScripts.RunBCalcScripts.calcBFromChunks import calcBFromChunks
from helperScripts.RunBCalcScripts.recmapHandler import calcRLengths
from helperScripts.RunBCalcScripts.recmapHandler import calcRDistances
import numpy as np

def process_single_chunk(chunk_num, chunk_size, blockstart, blockend, chr_start, chr_end,
                         calc_start, calc_end, num_chunks, precise_chunks,
                         lperchunk, b_values, rec_rate_per_chunk=None, silent=False):

    chunk_start = chr_start + chunk_num * chunk_size
    chunk_end   = min(chunk_start + chunk_size, calc_end)


    chunk_slice = b_values[chunk_start - calc_start:chunk_end - calc_start] # Get b_values for this chunk
    not_nan_mask = ~np.isnan(chunk_slice) # Make a mask for positions that are NOT NaN
    
    if not np.any(not_nan_mask):
        if not silent: print(f"No neutral sites in chunk {chunk_num}: {chunk_start}-{chunk_end}")
        return b_values
    
    pos_chunk = np.arange(chunk_start, chunk_end) # Array of positions in this chunk
    pos_chunk_clean   = pos_chunk[not_nan_mask] # Positions of neutral sites
    chunk_slice_clean = chunk_slice[not_nan_mask]

    B_from_distant_chunks = calcBFromChunks( # Compute B from distant chunks in non-precise region
        chunk_num, chunk_size,
        blockstart, blockend,
        chr_start, chr_end,
        calc_start, calc_end,
        num_chunks, precise_chunks,
        lperchunk,
        rec_rate_per_chunk
    )

    # Identify blocks in the "precise region"
    precise_region_start = np.maximum(chr_start, chr_start + (chunk_num - precise_chunks) * chunk_size)
    precise_region_end   = np.minimum(chr_end, chr_start + (chunk_num + 1 + precise_chunks) * chunk_size - 1)
    precise_blockregion_mask = (
        (precise_region_end   > blockstart) &
        (precise_region_start < blockend)
    )
    precise_blockstart = np.clip(blockstart[precise_blockregion_mask],
                                 a_min=precise_region_start, a_max=precise_region_end)
    precise_blockend   = np.clip(blockend[precise_blockregion_mask],
                                 a_min=precise_region_start, a_max=precise_region_end)

    physical_distances_upstream   = pos_chunk_clean[None, :] - precise_blockend[:, None] # All distances to blockends (upstream and downstream)
    physical_distances_downstream = precise_blockstart[:, None] - pos_chunk_clean[None, :] # All distances to blockstarts (upstream and downstream)

    downstream_mask = (pos_chunk_clean < precise_blockstart[:, None])

    upstream_mask   = (pos_chunk_clean > precise_blockend[:, None])
    flanking_mask   = downstream_mask | upstream_mask
    unique_indices, inverse_indices = np.unique(np.where(flanking_mask)[1], return_inverse=True)


    physical_distances = np.where( # Filter so only distances to blockends upstream and blockstarts downstream kept
        flanking_mask,
        np.where(upstream_mask, physical_distances_upstream, physical_distances_downstream),
        np.nan
    )
    flat_distances = physical_distances[flanking_mask] # Flatten array

    physical_lengths = precise_blockend - precise_blockstart
    flat_lengths   = np.repeat(physical_lengths, flanking_mask.sum(axis=1))
    nonzero_mask = flat_lengths != 0 # Remove genes of length 0
    flat_distances = flat_distances[nonzero_mask]
    flat_lengths   = flat_lengths[nonzero_mask]
    
    if rec_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE 
        precise_rates = rec_rate_per_chunk[np.maximum(0, chunk_num - precise_chunks):np.minimum(num_chunks, chunk_num + precise_chunks + 1)]
        rec_lengths = calcRLengths(precise_blockstart, precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, chunk_num)
        rec_distances_upstream, rec_distances_downstream = calcRDistances(precise_blockstart, precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, pos_chunk_clean, chunk_num, chunk_start)
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
    else:
        flank_B = calculateB_linear(flat_distances, flat_lengths)

    if flat_distances.size == 0 or flat_lengths.size == 0:
        flank_B = 1
        print(f"No nearby sites under selection in flanking region for chunk:", chunk_num)

    unique_indices, inverse_indices = np.unique(np.where(flanking_mask)[1], return_inverse=True)
    aggregated_B = np.ones_like(unique_indices, dtype=np.float64)
    np.multiply.at(aggregated_B, inverse_indices, flank_B) # Multiplicative sum of B calculated at a given site from multiple elements

    chunk_slice_clean[unique_indices] *= (aggregated_B * B_from_distant_chunks) # Update chunk slice and combine flank_B with B from distant chunks
    chunk_slice[not_nan_mask] = chunk_slice_clean # Put the updated (non-NaN) slice back into the original b_values
    mean_chunk_b = np.nanmean(chunk_slice) # Mean B for chunk

    if not silent: 
        print(f"Processing chunk: {pos_chunk.min()} - {pos_chunk.max()}")
        if rec_rate_per_chunk is not None:
            rec_rate = rec_rate_per_chunk[chunk_num]
            print(f"Chunk {chunk_num}: recombination rate = {rec_rate}")
        print(f"B from distant chunks: {B_from_distant_chunks}")
        print(f"Number of relevant genes: {len(precise_blockstart)}")
        print(f"Number of neutral sites in chunk [{chunk_start}-{chunk_end}): {np.isnan(chunk_slice).sum()}")
        print(f"Aggregated B values for chunk: {aggregated_B}")
        print(f"Mean B value for chunk {chunk_num}: [{chunk_start}-{chunk_end}]: {mean_chunk_b}")

    return b_values