import numpy as np
from helperScripts.calculate_B import calculate_B
from helperScripts.calcBFromChunks import calcBFromChunks
from multiprocessing import shared_memory

def process_single_chunk(chunk_num, chunk_size, blockstart, blockend, chr_start, chr_end, num_chunks, precise_chunks, lperchunk, b_values):

    # Generate positions for specified chunk
    chunk_start = chr_start + chunk_num * chunk_size
    chunk_end = min(chunk_start + chunk_size, chr_end)
    pos_chunk = np.arange(chunk_start, chunk_end)

    # Handle distant chunks for which l is combined
    precise_chunks = 3 # Controls how many adjacent chunks to the focal chunk have B calculated precisely rather than combined for all genes
    B_from_distant_chunks = calcBFromChunks(chunk_num, chunk_size, blockstart, blockend, chr_start, chr_end, num_chunks, precise_chunks, lperchunk)
    precise_region_start = chr_start + (chunk_num - precise_chunks) * chunk_size
    precise_region_end = chr_start + (chunk_num + precise_chunks + 1) * chunk_size

    # Find blockregions where both are within precise_region range
    precise_blockregion_mask = (precise_region_end >= blockstart) & (precise_region_start <= blockend)
    precise_blockstart = np.clip(blockstart[precise_blockregion_mask], a_min=precise_region_start, a_max=precise_region_end)
    precise_blockend = np.clip(blockend[precise_blockregion_mask], a_min=precise_region_start, a_max=precise_region_end)

    # Filter lengths to match the relevant blocks
    precise_lengths = precise_blockend - precise_blockstart

    # Calculate distances and masks for this chunk
    distances_downstream = precise_blockstart[:, None] - pos_chunk[None, :] 
    # print(distances_ups)
    distances_upstream = pos_chunk[None, :] - precise_blockend[:, None] 

    # Masks for flanking sites
    downstream_mask = (pos_chunk < precise_blockstart[:, None])
    upstream_mask = (pos_chunk > precise_blockend[:, None])
    flanking_mask = upstream_mask | downstream_mask

    # Combine the distances into a single array
    distances = np.where(flanking_mask, 
                        np.where(upstream_mask, distances_upstream, distances_downstream), 
                        np.nan)

    # Flatten distances and flanking_mask to match the selected elements
    flat_distances = distances[flanking_mask]  # Select distances where mask is True
    flat_lengths = np.repeat(precise_lengths, flanking_mask.sum(axis=1))[:len(flat_distances)]  # Repeat each length to match the mask, handling edge cases where lengths overshoot
    # Calculate B for the flattened data
    flank_B = calculate_B(flat_distances, flat_lengths)
    # Flatten flanking_mask to get indices of True values
    true_indices = np.where(flanking_mask)
    # Find unique indices and map each to its position in the unique list
    unique_indices, inverse_indices = np.unique(true_indices[1], return_inverse=True)
    # Aggregate B values for unique indices
    aggregated_B = np.ones_like(unique_indices, dtype=np.float64)
    np.multiply.at(aggregated_B, inverse_indices, flank_B)
    # Find global indices for the current chunk
    global_indices = pos_chunk[unique_indices] - chr_start  # Convert chunk indices to global indices
    # Vectorized updates for the `B` values in the b_values array
    b_values[global_indices] *= aggregated_B * B_from_distant_chunks


    print(f"Processing chunk: {pos_chunk.min()} - {pos_chunk.max()}")
    print(f"B from distant chunks: {B_from_distant_chunks}")
    print(f"Number of relevant genes: {len(precise_blockstart)}")
    print(f"Relevant blocks: {precise_blockstart}, {precise_blockend}")
    print(f"Aggregated B values for chunk: {aggregated_B}")

    return b_values