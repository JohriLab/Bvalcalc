import numpy as np

def extend_hri_regions_correction(b_values, rec_rate_per_chunk, hri_r_threshold = 0.1):

    ## First need to export properly from process_single_chunk.py

    print(f"Extending HRI regions", b_values, rec_rate_per_chunk)

    low_rec_chunk_ids = rec_rate_per_chunk < hri_r_threshold

    mask = low_rec_chunk_ids
    if mask.size == 0:
        return b_values

    # Inclusive starts/ends of each True-run
    interference_region_starts_idx = np.where(mask & np.r_[True, ~mask[:-1]])[0]
    interference_region_ends_idx   = np.where(mask & np.r_[~mask[1:], True])[0]

    print(interference_region_starts_idx, interference_region_ends_idx)

    return b_values