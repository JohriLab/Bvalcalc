from Bvalcalc.core.calculateB import calculateB_linear
import numpy as np

def calc_B_precise_noninterfering(precise_blockstart, precise_blockend, pos_chunk, chr_start, 
                                  chunk_size, chr_size, local_interference_indices):
    """
    Return B for the precise region *excluding* local interference blocks.
    """
    # 1) convert interfering chunk indices → bp range
    start_bp = chr_start + local_interference_indices[0] * chunk_size
    end_bp   = min(
        chr_start + (local_interference_indices[-1] + 1) * chunk_size - 1,
        chr_size
    )

    # 2) trim each block to drop overlap with [start_bp..end_bp]
    ts, te = [], []
    for bs, be in zip(precise_blockstart, precise_blockend):
        if be < start_bp or bs > end_bp:
            ts.append(bs);      te.append(be)       # fully outside
        else:
            if bs < start_bp:
                ts.append(bs);  te.append(start_bp - 1)  # left fragment
            if be > end_bp:
                ts.append(end_bp + 1); te.append(be)      # right fragment

    if not ts:
        return 1.0  # nothing remains → B = 1

    bs_arr = np.array(ts, dtype=int)
    be_arr = np.array(te, dtype=int)

    # 3) compute distance matrices
    up_dist   = pos_chunk[None, :] - be_arr[:, None]
    down_dist = bs_arr[:, None] - pos_chunk[None, :]
    up_mask   = pos_chunk >  be_arr[:, None]
    down_mask = pos_chunk <  bs_arr[:, None]
    mask      = up_mask | down_mask

    phys = np.where(mask, np.where(up_mask, up_dist, down_dist), np.nan)
    flat_distances = phys[mask]

    # 4) flatten lengths
    lengths = be_arr - bs_arr
    counts  = mask.sum(axis=1)
    flat_lengths = np.repeat(lengths, counts)

    # 5) drop any zero‐length entries
    valid = flat_lengths > 0
    flat_distances = flat_distances[valid]
    flat_lengths   = flat_lengths[valid]

    B_noninterfering_in_precise_region = calculateB_linear(flat_distances, flat_lengths)

    print(f"HAIRII", B_noninterfering_in_precise_region)

    # 6) linear B on non-interfering pieces
    return B_noninterfering_in_precise_region








    # if rec_rate_per_chunk is not None and gc_rate_per_chunk is not None: # IF REC_RATE MAP IS AVAILABLE and GC IS AVAILABLE
    #     B_noninterfering_in_precise_region = calculateB_recmap(distance_to_element=flat_distances, length_of_element=flat_lengths, rec_distances=flat_rec_distances, 
    #                                 rec_lengths=flat_rec_lengths, gc_distances=flat_gc_distances, gc_lengths=flat_gc_lengths)
    # elif rec_rate_per_chunk is not None and gc_rate_per_chunk is None: # IF REC_RATE MAP IS AVAILABLE and GC NOT AVAILABLE
    #     B_noninterfering_in_precise_region = calculateB_recmap(distance_to_element=flat_distances, length_of_element=flat_lengths, 
    #                                 rec_distances=flat_rec_distances, rec_lengths=flat_rec_lengths)
    # elif rec_rate_per_chunk is None and gc_rate_per_chunk is not None: # IF REC_RATE MAP NOT AVAILABLE and GC IS AVAILALBE
    #     B_noninterfering_in_precise_region = calculateB_recmap(distance_to_element=flat_distances, length_of_element=flat_lengths, 
    #                                 rec_distances=None, rec_lengths=None, gc_distances=flat_gc_distances, gc_lengths=flat_gc_lengths)
    # else: # NO MAPS AVAILABLE