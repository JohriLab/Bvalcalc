import numpy as np
from helperScripts.calculateB import calculateB
from helperScripts.RunBCalcScripts.calcBFromChunks import calcBFromChunks

def process_single_chunk(chunk_num, chunk_size, blockstart, blockend,
                         chr_start, chr_end, num_chunks, precise_chunks,
                         lperchunk, b_values):
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
        return b_values
    
    # Filter positions to skip pre-labeled NaNs
    pos_chunk_clean   = pos_chunk[not_nan_mask]
    chunk_slice_clean = chunk_slice[not_nan_mask]

    # == 1) If you rely on B_from_distant_chunks, compute that as usual ==
    B_from_distant_chunks = calcBFromChunks(
        chunk_num, chunk_size,
        blockstart, blockend,
        chr_start, chr_end,
        num_chunks, precise_chunks,
        lperchunk
    )

    # == 2) Identify blocks in the "precise region" (unchanged) ==
    precise_region_start = chr_start + (chunk_num - precise_chunks) * chunk_size
    precise_region_end   = chr_start + (chunk_num + precise_chunks + 1) * chunk_size

    precise_blockregion_mask = (
        (precise_region_end   >= blockstart) &
        (precise_region_start <= blockend)
    )
    precise_blockstart = np.clip(blockstart[precise_blockregion_mask],
                                 a_min=precise_region_start, a_max=precise_region_end)
    precise_blockend   = np.clip(blockend[precise_blockregion_mask],
                                 a_min=precise_region_start, a_max=precise_region_end)
    precise_lengths = precise_blockend - precise_blockstart

    # == 3) Do distance calculations ONLY for non-NaN sites ==
    distances_downstream = precise_blockstart[:, None] - pos_chunk_clean[None, :]
    distances_upstream   = pos_chunk_clean[None, :] - precise_blockend[:, None]

    downstream_mask = (pos_chunk_clean < precise_blockstart[:, None])
    upstream_mask   = (pos_chunk_clean > precise_blockend[:, None])
    flanking_mask   = downstream_mask | upstream_mask

    distances = np.where(
        flanking_mask,
        np.where(upstream_mask, distances_upstream, distances_downstream),
        np.nan
    )

    flat_distances = distances[flanking_mask]
    flat_lengths   = np.repeat(precise_lengths, flanking_mask.sum(axis=1))

    flank_B = calculateB(flat_distances, flat_lengths)

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

    print(f"Processing chunk: {pos_chunk.min()} - {pos_chunk.max()}")
    print(f"B from distant chunks: {B_from_distant_chunks}")
    print(f"Number of relevant genes: {len(precise_blockstart)}")
    print(f"Relevant blocks: {precise_blockstart}, {precise_blockend}")
    print(f"Number of NaN sites in chunk [{chunk_start}-{chunk_end}): {np.isnan(chunk_slice).sum()}")
    print(f"Aggregated B values for chunk: {aggregated_B}")
    print(f"Mean B value for chunk: [{chunk_start}-{chunk_end}]: {mean_chunk_b}")

    return b_values