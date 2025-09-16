import numpy as np
from Bvalcalc.core.helpers.calc_B_from_chunks import calc_B_from_chunks
from Bvalcalc.core.helpers.calc_B_precise_noninterfering import calc_B_precise_noninterfering
from Bvalcalc.core.calculateB import calculateB_hri

def calc_B_in_hri_region(quiet, chunk_idx, rec_rate_per_chunk, hri_r_threshold, lperchunk, chunk_size, chr_start, chr_size, num_chunks, gc_rate_per_chunk, precise_chunks, precise_blockstart, precise_blockend, pos_chunk, chunk_end, precise_region_start, precise_region_end, unlinked_B):
    if not quiet: print(f"Chunk {chunk_idx}: r modifier = {rec_rate_per_chunk[chunk_idx]}, which is at or below 0.1 threshold. Calculating B'. HRI calculation is enabled with --hri")
    low_rec_chunk_ids = rec_rate_per_chunk < hri_r_threshold

    # extend left, then right til there's a chunk with recombination
    interference_region_start_idx, interference_region_end_idx = chunk_idx, chunk_idx
    while interference_region_start_idx > 0 and low_rec_chunk_ids[interference_region_start_idx - 1]:
        interference_region_start_idx -= 1
    while interference_region_end_idx < low_rec_chunk_ids.size - 1 and low_rec_chunk_ids[interference_region_end_idx + 1]:
        interference_region_end_idx += 1
    total_interfering_L = lperchunk[interference_region_start_idx : interference_region_end_idx + 1].sum()## get combined 
    # if quiet: print(f"Contiguous low‐rec run: chunks {interference_region_start_idx}–{interference_region_end_idx}, interfering sites L = {total_interfering_L:.1f}")

    local_interference_indices = np.arange(interference_region_start_idx, interference_region_end_idx + 1) # Indexes for chunks in local interfering region

    B_from_distant_chunks = calc_B_from_chunks( # Re-compute B from distant chunks in non-precise region, exluding local interfering region
        chunk_idx, chunk_size, chr_start, chr_size, num_chunks, 
        precise_chunks, lperchunk, rec_rate_per_chunk, gc_rate_per_chunk, local_interference_indices)
    # print(f"Chunk {chunk_idx}: B_from_distant_chunks, excluding local interference region", B_from_distant_chunks)
    
    B_noninterfering_in_precise_region = calc_B_precise_noninterfering(precise_blockstart, precise_blockend, pos_chunk,
                                                                    chr_start, chunk_end, chunk_size, chr_size, precise_region_start, precise_region_end,local_interference_indices, chunk_idx, 
                                                                    rec_rate_per_chunk, gc_rate_per_chunk, quiet)

    U_lengths_in_low_rec_chunks = lperchunk[low_rec_chunk_ids]

    combined_prior_B = B_from_distant_chunks * B_noninterfering_in_precise_region[0] * unlinked_B

    interference_Bvals_per_chunk = calculateB_hri(
        distant_B=combined_prior_B,
        interfering_L=U_lengths_in_low_rec_chunks
    )

    return interference_Bvals_per_chunk