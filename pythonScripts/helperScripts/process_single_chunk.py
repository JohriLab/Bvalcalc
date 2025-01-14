import numpy as np
from helperScripts.calculate_B import calculate_B
from multiprocessing import shared_memory

def process_single_chunk(chunk_num, chunk_size, flank_blockstart, flank_blockend, blockstart, blockend, lengths, chr_start, chr_end, b_values):

    # Generate positions for specified chunk
    chunk_start = chr_start + chunk_num * chunk_size
    chunk_end = min(chunk_start + chunk_size, chr_end)
    pos_chunk = np.arange(chunk_start, chunk_end)

    # Filter blockstart and blockend to only include where flanking region overlaps with chunk
    relevant_blockregion = (pos_chunk.max() >= flank_blockstart) & (pos_chunk.min() <= flank_blockend)
    relevant_blockstart = blockstart[relevant_blockregion]
    relevant_blockend = blockend[relevant_blockregion]
    relevant_flank_blockstart = flank_blockstart[relevant_blockregion]
    relevant_flank_blockend = flank_blockend[relevant_blockregion]

    # Filter lengths to match the relevant blocks
    relevant_lengths = lengths[relevant_blockregion]

    # Calculate distances and masks for this chunk
    distances_upstream = relevant_blockstart[:, None] - pos_chunk[None, :] 
    distances_downstream = pos_chunk[None, :] - relevant_blockend[:, None] 

    # Masks for flanking sites
    upstream_mask = (pos_chunk < relevant_blockstart[:, None]) & \
                    (pos_chunk > (relevant_flank_blockstart[:, None]))
    downstream_mask = (pos_chunk > relevant_blockend[:, None]) & \
                    (pos_chunk < (relevant_flank_blockend[:, None]))
    flanking_mask = upstream_mask | downstream_mask

    # Combine the distances into a single array
    distances = np.where(flanking_mask, 
                        np.where(upstream_mask, distances_upstream, distances_downstream), 
                        np.nan)
    # Flatten distances and flanking_mask to match the selected elements
    flat_distances = distances[flanking_mask]  # Select distances where mask is True
    flat_lengths = np.repeat(relevant_lengths, flanking_mask.sum(axis=1))[:len(flat_distances)]  # Repeat each length to match the mask, handling edge cases where lengths overshoot
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
    b_values[global_indices] *= aggregated_B

    print(f"Processing chunk: {pos_chunk.min()} - {pos_chunk.max()}")
    print(f"Number of relevant genes: {len(relevant_blockstart)}")
    print(f"Relevant blocks: {relevant_blockstart}, {relevant_blockend}")
    print(f"Aggregated B values for chunk: {aggregated_B}")

    return b_values