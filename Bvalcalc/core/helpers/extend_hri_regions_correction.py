import numpy as np

def extend_hri_regions_correction(b_values, rec_rate_per_chunk, chunk_size, chr_start, calc_start, calc_end, hri_r_threshold = 0.1):

    ## First need to export properly from process_single_chunk.py

    # print(f"Extending HRI regions", b_values, rec_rate_per_chunk)

    low_rec_chunk_ids = rec_rate_per_chunk < hri_r_threshold

    mask = low_rec_chunk_ids
    if not mask.any():
        return b_values

    # Inclusive starts/ends of each True-run (contiguous HRI regions)
    interference_region_starts_idx = np.where(mask & np.r_[True, ~mask[:-1]])[0]
    interference_region_ends_idx   = np.where(mask & np.r_[~mask[1:], True])[0]

    # print(interference_region_starts_idx, interference_region_ends_idx)

    base_chunk_idx = (calc_start - chr_start) // chunk_size

    # Get start and end positions of interference regions
    abs_start_chunks = base_chunk_idx + interference_region_starts_idx
    abs_end_chunks   = base_chunk_idx + interference_region_ends_idx + 1  # +1 for end-exclusive

    interference_region_start_pos = chr_start + abs_start_chunks * chunk_size
    interference_region_end_pos = np.minimum(calc_end,
        chr_start + abs_end_chunks * chunk_size - 1)

    # Filter to only include interference regions that overlap with the calculation region
    valid_mask = (interference_region_end_pos >= calc_start) & (interference_region_start_pos <= calc_end)
    if not valid_mask.any():
        return b_values

    interference_region_start_pos = interference_region_start_pos[valid_mask]
    interference_region_end_pos = interference_region_end_pos[valid_mask]

    # Clip positions to be within the calculation region
    interference_region_start_pos = np.maximum(interference_region_start_pos, calc_start)
    interference_region_end_pos = np.minimum(interference_region_end_pos, calc_end)

    B_in_interference_regions = b_values[interference_region_start_pos - calc_start]

    # Debug peek
    # print("interference_region_start_pos:", interference_region_start_pos)
    # print("interference_region_end_pos  :", interference_region_end_pos)
    # print("b_values  :", B_in_interference_regions)

    left_extended  = np.zeros(len(interference_region_start_pos), dtype=int)
    right_extended = np.zeros(len(interference_region_end_pos), dtype=int)

    n = b_values.shape[0]


    # ---- Backward extension from each region start ----
    # For region i, extend left from start_pos[i] - 1 while region B is higher
    # than the current b_values. Stop at calc_start or just after the previous region.
    for i in range(len(interference_region_start_pos)):
        b_inside = B_in_interference_regions[i]

        # ---- backward extension ----
        idx_left = (interference_region_start_pos[i] - calc_start) - 1
        if i == 0:
            stop_left_rel = 0
        else:
            stop_left_abs = interference_region_end_pos[i - 1] + 1
            stop_left_rel = max(0, stop_left_abs - calc_start)

        while idx_left >= stop_left_rel:
            if b_inside > b_values[idx_left]:
                b_values[idx_left] = b_inside
                left_extended[i] += 1
                # NEW: stop if extended beyond chunk_size
                if left_extended[i] > chunk_size:
                    break
                idx_left -= 1
            else:
                break

        # ---- forward extension ----
        idx_right = (interference_region_end_pos[i] - calc_start) + 1
        if i == len(interference_region_start_pos) - 1:
            stop_right_rel = n - 1
        else:
            stop_right_abs = interference_region_start_pos[i + 1] - 1
            stop_right_rel = min(n - 1, stop_right_abs - calc_start)

        while idx_right <= stop_right_rel:
            if b_inside > b_values[idx_right]:
                b_values[idx_right] = b_inside
                right_extended[i] += 1
                # NEW: stop if extended beyond chunk_size
                if right_extended[i] > chunk_size:
                    break
                idx_right += 1
            else:
                break

    # ---- warnings ----
    for i in range(len(interference_region_start_pos)):
        warn_left  = left_extended[i]  > chunk_size
        warn_right = right_extended[i] > chunk_size
        if warn_left or warn_right:
            sides = []
            if warn_left:
                sides.append(f"left={left_extended[i]}")
            if warn_right:
                sides.append(f"right={right_extended[i]}")
            print(
                f"WARNING: HRI extension exceeded chunk_size at region {i} "
                f"[{interference_region_start_pos[i]}-{interference_region_end_pos[i]}]; "
                f"{', '.join(sides)} bases > chunk_size ({chunk_size}). "
                f"This suggests B is very low and HRI is too pervasive for B to be calculated effectively in this region."
            )

        # NEW: compute final (start,end) positions after extension, inclusive
    extended_starts = np.maximum(interference_region_start_pos - left_extended, calc_start)
    extended_ends   = np.minimum(interference_region_end_pos + right_extended, calc_end)
    # extended_regions_pos = np.stack([extended_starts, extended_ends], axis=1).astype(int)

    print(f"Extended HRI regions to {extended_starts, extended_ends}")

    return extended_starts, extended_ends