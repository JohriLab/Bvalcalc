import numpy as np

def extend_hri_regions_correction(b_values, rec_rate_per_chunk, chunk_size, chr_start, calc_start, calc_end, hri_r_threshold = 0.1):

    ## First need to export properly from process_single_chunk.py

    print(f"Extending HRI regions", b_values, rec_rate_per_chunk)

    low_rec_chunk_ids = rec_rate_per_chunk < hri_r_threshold

    mask = low_rec_chunk_ids
    if not mask.any():
        return b_values

    # Inclusive starts/ends of each True-run (contiguous HRI regions)
    interference_region_starts_idx = np.where(mask & np.r_[True, ~mask[:-1]])[0]
    interference_region_ends_idx   = np.where(mask & np.r_[~mask[1:], True])[0]

    print(interference_region_starts_idx, interference_region_ends_idx)

    base_chunk_idx = (calc_start - chr_start) // chunk_size

    # Get start and end positions of interference regions
    abs_start_chunks = base_chunk_idx + interference_region_starts_idx
    abs_end_chunks   = base_chunk_idx + interference_region_ends_idx + 1  # +1 for end-exclusive

    interference_region_start_pos = chr_start + abs_start_chunks * chunk_size

    interference_region_end_pos = np.minimum(calc_end,
        chr_start + abs_end_chunks * chunk_size - 1)
    
    B_in_interference_regions = b_values[interference_region_start_pos-calc_start]

    # Debug peek
    print("interference_region_start_pos:", interference_region_start_pos)
    print("interference_region_end_pos  :", interference_region_end_pos)
    print("b_values  :", b_values[interference_region_start_pos-calc_start])


    return b_values