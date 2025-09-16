from Bvalcalc.core.calculateB import calculateB_linear, calculateB_recmap
from Bvalcalc.core.helpers.calc_B_from_chunks import calc_B_from_chunks
from Bvalcalc.core.helpers.calc_R_len_dist import calc_R_lengths
from Bvalcalc.core.helpers.calc_R_len_dist import calc_R_distances
from Bvalcalc.core.helpers.calc_B_in_genes import calc_B_in_genes
from Bvalcalc.core.helpers.calc_B_precise_noninterfering import calc_B_precise_noninterfering
import numpy as np

def process_single_chunk(chunk_idx, chunk_size, blockstart, blockend, chr_start, chr_size,
                         calc_start, calc_end, num_chunks, precise_chunks,lperchunk, 
                         b_values, rec_rate_per_chunk=None, gc_rate_per_chunk=None, hri=False, quiet=False, verbose=False, unlinked_B=1.0):
    
    chunk_start =  chr_start + chunk_idx * chunk_size
    chunk_end   = min(chunk_start + chunk_size - 1, calc_end)
    chunk_slice = b_values[chunk_start:chunk_end+1] # Get b_values for this chunk
    pos_chunk = np.arange(chunk_start, chunk_end+1) # Array of positions in this chunk
    
    # Identify blocks in the "precise region"
    precise_region_start = np.maximum(chr_start, chr_start + (chunk_idx - precise_chunks) * chunk_size)
    precise_region_end   = np.minimum(chr_size, chr_start + (chunk_idx + 1 + precise_chunks) * chunk_size - 1)
    precise_blockregion_mask = (
        (precise_region_end   > blockstart) &
        (precise_region_start < blockend))
    precise_blockstart = np.clip(blockstart[precise_blockregion_mask],a_min=precise_region_start, a_max=precise_region_end)
    precise_blockend   = np.clip(blockend[precise_blockregion_mask], a_min=precise_region_start, a_max=precise_region_end)

    physical_distances_upstream   = pos_chunk[None, :] - precise_blockend[:, None] # All distances to blockends (upstream and downstream)
    physical_distances_downstream = precise_blockstart[:, None] - pos_chunk[None, :] # All distances to blockstarts (upstream and downstream)

    downstream_mask = (pos_chunk < precise_blockstart[:, None]) # True when position is less than blockstart (gene is downstream) 
    upstream_mask   = (pos_chunk > precise_blockend[:, None]) # True when position is more than blockend (gene is upstream)
    flanking_mask   = downstream_mask | upstream_mask 

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
        precise_rates = rec_rate_per_chunk[np.maximum(0, chunk_idx - precise_chunks):np.minimum(num_chunks, chunk_idx + precise_chunks + 1)]
        rec_lengths = calc_R_lengths(precise_blockstart , precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size)
        rec_distances_upstream, rec_distances_downstream = calc_R_distances(precise_blockstart, precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, pos_chunk, chunk_start)
        rec_distances = np.where(
            flanking_mask,
            np.where(upstream_mask, rec_distances_upstream, rec_distances_downstream),
            np.nan)
        flat_rec_distances = rec_distances[flanking_mask]
        flat_rec_lengths   = np.repeat(rec_lengths, flanking_mask.sum(axis=1))
        nonzero_rec_mask = flat_lengths != 0 # Remove genes of length 0
        flat_rec_distances = flat_rec_distances[nonzero_rec_mask]
        flat_rec_lengths   = flat_rec_lengths[nonzero_rec_mask]

    if gc_rate_per_chunk is not None: # IF GC_RATE MAP IS AVAILABLE 
        precise_gc_rates = gc_rate_per_chunk[np.maximum(0, chunk_idx - precise_chunks):np.minimum(num_chunks, chunk_idx + precise_chunks + 1)]
        gc_lengths = calc_R_lengths(precise_blockstart, precise_blockend, precise_gc_rates, precise_region_start, precise_region_end, chunk_size)
        gc_distances_upstream, gc_distances_downstream = calc_R_distances(precise_blockstart, precise_blockend, precise_gc_rates, precise_region_start, precise_region_end, chunk_size, pos_chunk, chunk_start)
        gc_distances = np.where(
            flanking_mask,
            np.where(upstream_mask, gc_distances_upstream, gc_distances_downstream),
            np.nan)
        flat_gc_distances = gc_distances[flanking_mask]
        flat_gc_lengths   = np.repeat(gc_lengths, flanking_mask.sum(axis=1))
        nonzero_gc_mask = flat_lengths != 0 # Remove genes of length 0
        flat_gc_distances = flat_gc_distances[nonzero_gc_mask]
        flat_gc_lengths   = flat_gc_lengths[nonzero_gc_mask]

    B_from_distant_chunks = calc_B_from_chunks( # Compute B from distant chunks in non-precise region
        chunk_idx, chunk_size, chr_start, chr_size, num_chunks, 
        precise_chunks, lperchunk, rec_rate_per_chunk, gc_rate_per_chunk)
    
    # Calculate B for genes within chunk. within_gene_B is to include B for genic sites from BGS caused by the gene they're in
    if rec_rate_per_chunk is not None and gc_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE and GC IS AVAILABLE
        within_gene_B = calc_B_in_genes(chunk_size, num_chunks, precise_chunks, precise_blockstart, precise_blockend, chunk_start, chunk_end+1, physical_lengths, precise_region_start, chunk_idx, rec_rate_per_chunk = rec_rate_per_chunk, gc_rate_per_chunk = gc_rate_per_chunk, rec_lengths = rec_lengths, gc_lengths = gc_lengths)
        flank_B = calculateB_recmap(distance_to_element=flat_distances, length_of_element=flat_lengths, rec_distances=flat_rec_distances, 
                                    rec_lengths=flat_rec_lengths, gc_distances=flat_gc_distances, gc_lengths=flat_gc_lengths)
    elif rec_rate_per_chunk is not None and gc_rate_per_chunk is None: # IF REC_RATE MAP IS AVAILABLE and GC NOT AVAILABLE
        within_gene_B = calc_B_in_genes(chunk_size, num_chunks, precise_chunks, precise_blockstart, precise_blockend, chunk_start, chunk_end+1, physical_lengths, precise_region_start, chunk_idx, rec_rate_per_chunk = rec_rate_per_chunk, gc_rate_per_chunk = None, rec_lengths = rec_lengths, gc_lengths = None)
        flank_B = calculateB_recmap(distance_to_element=flat_distances, length_of_element=flat_lengths, 
                                    rec_distances=flat_rec_distances, rec_lengths=flat_rec_lengths)
    elif rec_rate_per_chunk is None and gc_rate_per_chunk is not None: # IF REC_RATE MAP NOT AVAILABLE and GC IS AVAILALBE
        within_gene_B = calc_B_in_genes(chunk_size, num_chunks, precise_chunks, precise_blockstart, precise_blockend, chunk_start, chunk_end+1, physical_lengths, precise_region_start, chunk_idx, rec_rate_per_chunk = None, gc_rate_per_chunk = gc_rate_per_chunk, rec_lengths = None, gc_lengths = gc_lengths)
        flank_B = calculateB_recmap(distance_to_element=flat_distances, length_of_element=flat_lengths, 
                                    rec_distances=None, rec_lengths=None, gc_distances=flat_gc_distances, gc_lengths=flat_gc_lengths)
    else: # NO MAPS AVAILABLE
        within_gene_B = calc_B_in_genes(chunk_size, num_chunks, precise_chunks, precise_blockstart, precise_blockend, chunk_start, chunk_end+1, physical_lengths, precise_region_start, chunk_idx, rec_rate_per_chunk = None, gc_rate_per_chunk = None, rec_lengths = None, gc_lengths = None)
        flank_B = calculateB_linear(flat_distances, flat_lengths)

    # Combine B's calculated from distant genes, and genes within the region!
    safe_flank_B = np.concatenate((np.ones(chunk_end+1 - chunk_start, dtype=float), flank_B)) # Add an array of flank_B where all sites are B = 1, to account for sites with no flanking genes
    new_flanking_mask = np.concatenate((np.ones((1, chunk_end+1 - chunk_start), dtype=bool), flanking_mask), axis=0)
    unique_indices, inverse_indices = np.unique(np.where(new_flanking_mask)[1], return_inverse=True)
    aggregated_B = np.ones_like(np.ones_like(np.arange(chunk_start,chunk_end+1), dtype=np.float64), dtype=np.float64)
    np.multiply.at(aggregated_B, inverse_indices, safe_flank_B) # Multiplicative sum of B calculated at a given site from multiple elements
    hri_r_threshold = 0.1 # fraction of "r" in a chunk that triggers Bprime hri calculation
    hri_L_threshold = 1000 # minimum number of selected sites in a chunk that triggers Bprime hri calculation
    
    # Check if this chunk should trigger HRI calculation
    should_do_hri = False
    if (hri and rec_rate_per_chunk is not None and rec_rate_per_chunk[chunk_idx] < hri_r_threshold):
        # Calculate the total L across the entire interference region
        low_rec_chunk_ids = rec_rate_per_chunk < hri_r_threshold
        interference_region_start_idx, interference_region_end_idx = chunk_idx, chunk_idx
        while interference_region_start_idx > 0 and low_rec_chunk_ids[interference_region_start_idx - 1]:
            interference_region_start_idx -= 1
        while interference_region_end_idx < low_rec_chunk_ids.size - 1 and low_rec_chunk_ids[interference_region_end_idx + 1]:
            interference_region_end_idx += 1
        total_interfering_L = lperchunk[interference_region_start_idx : interference_region_end_idx + 1].sum()
        should_do_hri = total_interfering_L > hri_L_threshold
    
    if should_do_hri: # Only do this if user has --hri active
        from Bvalcalc.core.helpers.calc_B_in_hri_region import calc_B_in_hri_region
        hri_aggregated_B = calc_B_in_hri_region(quiet, chunk_idx, rec_rate_per_chunk, hri_r_threshold, lperchunk, chunk_size, chr_start, chr_size, num_chunks, gc_rate_per_chunk, precise_chunks, precise_blockstart, precise_blockend, pos_chunk, chunk_end, precise_region_start, precise_region_end, unlinked_B)
        chunk_slice *= hri_aggregated_B
        return b_values
    else:
        if unique_indices.size == 0: # If there are no nearby sites under selection
            chunk_slice *= (B_from_distant_chunks * within_gene_B)
        else:
            chunk_slice *= (aggregated_B * B_from_distant_chunks * within_gene_B) # Update chunk slice and combine flank_B with B from distant chunks

        mean_chunk_b = np.nanmean(chunk_slice) # Mean B for chunk




    if verbose: # Per-chunk summaries
        print(f"Processing chunk {chunk_idx}: {pos_chunk.min()} - {pos_chunk.max()}")
        if rec_rate_per_chunk is not None:
            print(f"Chunk {chunk_idx}: recombination rate = {rec_rate_per_chunk[chunk_idx]}")
        if gc_rate_per_chunk is not None:
            print(f"Chunk {chunk_idx}: gene conversion rate = {gc_rate_per_chunk[chunk_idx]}")
        # print(f"B from distant chunks: {B_from_distant_chunks}")
        print(f"Number of relevant genes: {len(precise_blockstart)}")
        print(f"Number of neutral sites in chunk [{chunk_start}-{chunk_end}): {np.isnan(chunk_slice).sum()}")
        # print(f"Aggregated B values for chunk: {aggregated_B}")
        print(f"Mean B value for chunk {chunk_idx}: [{chunk_start}-{chunk_end}]: {mean_chunk_b}")

    return b_values