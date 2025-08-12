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
    
    B_in_interference_regions = b_values[np.maximum(interference_region_start_pos-calc_start, 0)] # Find in b_values array which starts at calc_start

    # Debug peek
    print("interference_region_start_pos:", interference_region_start_pos)
    print("interference_region_end_pos  :", interference_region_end_pos)
    print("b_values  :", B_in_interference_regions)


    # ---- Backward extension from each region start ----
    # For region i, extend left from start_pos[i] - 1 while region B is higher
    # than the current b_values. Stop at calc_start or just after the previous region.
    for i in range(len(interference_region_start_pos)):
        b_inside = B_in_interference_regions[i]

        # index in b_values for the neighbor left of the region start
        idx = (interference_region_start_pos[i] - calc_start) - 1
        if idx < 0:
            continue  # already at the leftmost calc position

        # Do not overwrite into the previous interference region
        if i == 0:
            left_stop_rel = 0  # calc_start in b_values coordinates
        else:
            # previous region's end is inclusive; stop just to its right
            left_stop_abs = interference_region_end_pos[i - 1] + 1
            left_stop_rel = max(0, left_stop_abs - calc_start)

        # Walk left while region B is higher than current b_values
        while idx >= left_stop_rel:
            if b_inside > b_values[idx]:
                b_values[idx] = b_inside
                idx -= 1
            else:
                break

    # ---- Forward extension from each region end ----
    # Mirror of above: extend right from end_pos[i] + 1 while region B is higher.
    # Stop at calc_end or just before the next region's start.
    n = b_values.shape[0]
    for i in range(len(interference_region_end_pos)):
        b_inside = B_in_interference_regions[i]

        # index in b_values for the neighbor right of the region end
        idx = (interference_region_end_pos[i] - calc_start) + 1
        if idx >= n:
            continue  # already at or beyond the rightmost calc position

        # Do not overwrite into the next interference region
        if i == len(interference_region_start_pos) - 1:
            right_stop_rel = n - 1  # last valid index
        else:
            # next region's start is absolute; stop just to its left
            right_stop_abs = interference_region_start_pos[i + 1] - 1
            right_stop_rel = min(n - 1, right_stop_abs - calc_start)

        # Walk right while region B is higher than current b_values
        while idx <= right_stop_rel:
            if b_inside > b_values[idx]:
                b_values[idx] = b_inside
                idx += 1
            else:
                break

    return b_values