from Bvalcalc.core.calculateB import calculateB_linear, calculateB_recmap
from Bvalcalc.core.helpers.calc_R_len_dist import calc_R_lengths, calc_R_distances
import numpy as np

def calc_B_precise_noninterfering(
    precise_blockstart, precise_blockend, pos_chunk,
    chr_start, chunk_size, chr_size,
    local_interference_indices, chunk_idx,
    rec_rate_per_chunk, gc_rate_per_chunk
):
    """
    Calculate B in the precise region excluding local interference blocks,
    with optional rec/GCall formulas.

    DEV note: Lots of code is redundant here with the parent process_single_chunk(). 
    DEV note: Could be streamlined but given HRI calculation is assumed to be relatively rare, the speed hit should be manageable.
    """
    # 1) bp range of interference region
    start_bp = chr_start + local_interference_indices[0] * chunk_size
    end_bp   = min(
        chr_start + (local_interference_indices[-1] + 1) * chunk_size - 1,
        chr_size
    )
    print(f"[DEBUG] interference bp range: {start_bp}-{end_bp}")

    # 2) trim blocks to remove overlap
    ts, te = [], []
    for bs, be in zip(precise_blockstart, precise_blockend):
        if be < start_bp or bs > end_bp:
            ts.append(bs);      te.append(be)       # keep whole
        else:
            if bs < start_bp:
                ts.append(bs);  te.append(start_bp - 1)  # left fragment
            if be > end_bp:
                ts.append(end_bp + 1); te.append(be)      # right fragment
    print(f"[DEBUG] original blocks: {list(zip(precise_blockstart, precise_blockend))}")
    if not ts:
        print(f"[DEBUG] no non‐interfering blocks → B=1")
        return 1.0
    bs_arr = np.array(ts, dtype=int)
    be_arr = np.array(te, dtype=int)
    print(f"[DEBUG] trimmed blocks:  {list(zip(bs_arr, be_arr))}")

    # 3) build distance mask
    up_dist   = pos_chunk[None, :] - be_arr[:, None]
    down_dist = bs_arr[:, None] - pos_chunk[None, :]
    up_mask   = pos_chunk >  be_arr[:, None]
    down_mask = pos_chunk <  bs_arr[:, None]
    mask      = up_mask | down_mask

    phys = np.where(mask, np.where(up_mask, up_dist, down_dist), np.nan)
    flat_distances = phys[mask]
    print(f"[DEBUG] flat_distances ({len(flat_distances)}): {flat_distances[:5]}{'...' if len(flat_distances)>5 else ''}")

    # 4) flatten lengths
    lengths = be_arr - bs_arr
    counts  = mask.sum(axis=1)
    flat_lengths = np.repeat(lengths, counts)
    print(f"[DEBUG] flat_lengths   ({len(flat_lengths)}): {flat_lengths[:5]}{'...' if len(flat_lengths)>5 else ''}")

    # 5) drop zeros
    valid = flat_lengths > 0
    flat_distances = flat_distances[valid]
    flat_lengths   = flat_lengths[valid]

    # 6) optionally flatten rec
    if rec_rate_per_chunk is not None:
        # slice recomb rates to precise region
        region_chunk_start = (start_bp - chr_start) // chunk_size
        region_chunk_end   = (end_bp   - chr_start) // chunk_size
        r_rates = rec_rate_per_chunk[region_chunk_start:region_chunk_end+1]

        # weighted block lengths
        rec_lens = calc_R_lengths(
            bs_arr, be_arr,
            r_rates,
            start_bp, end_bp,
            chunk_size
        )

        # weighted distances
        focal_start_bp = chr_start + chunk_idx * chunk_size
        rec_up, rec_down = calc_R_distances(
            bs_arr, be_arr,
            r_rates,
            start_bp, end_bp,
            chunk_size,
            pos_chunk,
            focal_start_bp
        )
        rec_phys = np.where(mask, np.where(up_mask, rec_up, rec_down), np.nan)
        flat_rec_distances = rec_phys[mask][valid]
        flat_rec_lengths   = np.repeat(rec_lens, counts)[valid]
        print(f"[DEBUG] flat_rec_distances ({len(flat_rec_distances)}): {flat_rec_distances[:5]}{'...' if len(flat_rec_distances)>5 else ''}")
        print(f"[DEBUG] flat_rec_lengths   ({len(flat_rec_lengths)}): {flat_rec_lengths[:5]}{'...' if len(flat_rec_lengths)>5 else ''}")
    else:
        flat_rec_distances = flat_rec_lengths = None

    # 7) optionally flatten gc
    if gc_rate_per_chunk is not None:
        # slice gc rates to precise region
        region_chunk_start = (start_bp - chr_start) // chunk_size
        region_chunk_end   = (end_bp   - chr_start) // chunk_size
        g_rates = gc_rate_per_chunk[region_chunk_start:region_chunk_end+1]

        # weighted block lengths
        gc_lens = calc_R_lengths(
            bs_arr, be_arr,
            g_rates,
            start_bp, end_bp,
            chunk_size
        )

        # weighted distances
        focal_start_bp = chr_start + chunk_idx * chunk_size
        gc_up, gc_down = calc_R_distances(
            bs_arr, be_arr,
            g_rates,
            start_bp, end_bp,
            chunk_size,
            pos_chunk,
            focal_start_bp
        )
        gc_phys = np.where(mask, np.where(up_mask, gc_up, gc_down), np.nan)
        flat_gc_distances = gc_phys[mask][valid]
        flat_gc_lengths   = np.repeat(gc_lens, counts)[valid]
        print(f"[DEBUG] flat_gc_distances ({len(flat_gc_distances)}): {flat_gc_distances[:5]}{'...' if len(flat_gc_distances)>5 else ''}")
        print(f"[DEBUG] flat_gc_lengths   ({len(flat_gc_lengths)}): {flat_gc_lengths[:5]}{'...' if len(flat_gc_lengths)>5 else ''}")
    else:
        flat_gc_distances = flat_gc_lengths = None

    # 8) dispatch to the correct B function
    if flat_rec_distances is not None and flat_gc_distances is not None:
        B = calculateB_recmap(
            distance_to_element=flat_distances,
            length_of_element=flat_lengths,
            rec_distances=flat_rec_distances,
            rec_lengths=flat_rec_lengths,
            gc_distances=flat_gc_distances,
            gc_lengths=flat_gc_lengths
        )
    elif flat_rec_distances is not None:
        B = calculateB_recmap(
            distance_to_element=flat_distances,
            length_of_element=flat_lengths,
            rec_distances=flat_rec_distances,
            rec_lengths=flat_rec_lengths
        )
    elif flat_gc_distances is not None:
        B = calculateB_recmap(
            distance_to_element=flat_distances,
            length_of_element=flat_lengths,
            gc_distances=flat_gc_distances,
            gc_lengths=flat_gc_lengths
        )
    else:
        B = calculateB_linear(flat_distances, flat_lengths)

    print(f"[DEBUG] final B_noninterfering = {B}")
    return B
