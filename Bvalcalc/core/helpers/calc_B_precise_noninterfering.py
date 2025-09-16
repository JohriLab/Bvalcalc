from Bvalcalc.core.calculateB import calculateB_linear, calculateB_recmap
from Bvalcalc.core.helpers.calc_R_len_dist import calc_R_lengths, calc_R_distances
import numpy as np

def calc_B_precise_noninterfering(
    precise_blockstart, precise_blockend, pos_chunk,
    chr_start, chunk_end, chunk_size, chr_size, precise_region_start, precise_region_end,
    local_interference_indices, chunk_idx,
    rec_rate_per_chunk, gc_rate_per_chunk, quiet=False
):
    """
    Calculate B in the precise region excluding local interference blocks,
    with optional recombination/gene-conversion formulas.

    DEV note: Lots of code is redundant here with the parent process_single_chunk().
    DEV note: Could be streamlined but given HRI calculation is assumed to be relatively
    rare, the speed hit should be manageable.
    """
    # 1) bp range of interference region
    start_bp = chr_start + local_interference_indices[0] * chunk_size
    end_bp   = min(
        chr_start + (local_interference_indices[-1] + 1) * chunk_size - 1,
        chr_size
    )
    # print(f"[DEBUG] interference bp range: {start_bp}-{end_bp}")

    # 2) trim blocks to remove overlap
    ts, te = [], []
    for bs, be in zip(precise_blockstart, precise_blockend):
        if be < start_bp or bs > end_bp:
            ts.append(bs); te.append(be)
        else:
            if bs < start_bp:
                ts.append(bs); te.append(start_bp - 1)
            if be > end_bp:
                ts.append(end_bp + 1); te.append(be)
    # print(f"[DEBUG] original blocks: {list(zip(precise_blockstart, precise_blockend))}")
    # if not ts:
    #     print(f"Chunk {chunk_idx}: No non-interfering blocks; B=1")
    #     return 1.0
    bs_arr = np.array(ts, dtype=int)
    be_arr = np.array(te, dtype=int)
    # print(f"[DEBUG] trimmed blocks:  {list(zip(bs_arr, be_arr))}")

    # 3) build distance mask
    up_dist   = pos_chunk[None, :] - be_arr[:, None]
    down_dist = bs_arr[:, None] - pos_chunk[None, :]
    up_mask   = pos_chunk > be_arr[:, None]
    down_mask = pos_chunk < bs_arr[:, None]
    mask      = up_mask | down_mask

    phys = np.where(mask, np.where(up_mask, up_dist, down_dist), np.nan)
    flat_distances = phys[mask]

    # print(f"[DEBUG] flat_distances ({len(flat_distances)}): "
    #       f"{flat_distances[:5]}{'...' if len(flat_distances)>5 else ''}")

    # 4) flatten lengths
    lengths = be_arr - bs_arr
    counts  = mask.sum(axis=1)
    flat_lengths = np.repeat(lengths, counts)
    # print(f"[DEBUG] flat_lengths   ({len(flat_lengths)}): "
    #       f"{flat_lengths[:5]}{'...' if len(flat_lengths)>5 else ''}")

    # 5) drop zeros
    valid = flat_lengths > 0
    flat_distances = flat_distances[valid]
    flat_lengths   = flat_lengths[valid]

    # 6) optionally flatten recombination
    if rec_rate_per_chunk is not None:
        region_start_idx = (precise_region_start.min() - chr_start) // chunk_size
        region_end_idx   = (precise_region_end.max() - chr_start) // chunk_size
        r_rates = rec_rate_per_chunk[region_start_idx:region_end_idx+1]

        rec_lens = calc_R_lengths(
            bs_arr, be_arr,
            r_rates,
            precise_region_start, precise_region_end,
            chunk_size
        )
        chunk_start = chr_start + chunk_idx * chunk_size
        rec_up, rec_down = calc_R_distances(
            bs_arr, be_arr,
            r_rates,
            precise_region_start, precise_region_end,
            chunk_size,
            pos_chunk,
            chunk_start
        )
        rec_phys = np.where(mask, np.where(up_mask, rec_up, rec_down), np.nan)
        flat_rec_distances = rec_phys[mask][valid]
        flat_rec_lengths   = np.repeat(rec_lens, counts)[valid]
        # print(f"[DEBUG] flat_rec_distances ({len(flat_rec_distances)}): "
        #       f"{flat_rec_distances[:5]}{'...' if len(flat_rec_distances)>5 else ''}")
        # print(f"[DEBUG] flat_rec_lengths    ({len(flat_rec_lengths)}):  "
        #       f"{flat_rec_lengths[:5]}{'...' if len(flat_rec_lengths)>5 else ''}")
    else:
        flat_rec_distances = flat_rec_lengths = None

    # 7) optionally flatten gene-conversion
    if gc_rate_per_chunk is not None:
        region_start_idx = (precise_region_start - chr_start) // chunk_size
        region_end_idx   = (precise_region_end - chr_start) // chunk_size
        g_rates = gc_rate_per_chunk[region_start_idx:region_end_idx+1]

        gc_lens = calc_R_lengths(
            bs_arr, be_arr,
            g_rates,
            precise_region_start, precise_region_end,
            chunk_size
        )
        chunk_start = chr_start + chunk_idx * chunk_size
        gc_up, gc_down = calc_R_distances(
            bs_arr, be_arr,
            g_rates,
            precise_region_start, precise_region_end,
            chunk_size,
            pos_chunk,
            chunk_start
        )
        gc_phys = np.where(mask, np.where(up_mask, gc_up, gc_down), np.nan)
        flat_gc_distances = gc_phys[mask][valid]
        flat_gc_lengths   = np.repeat(gc_lens, counts)[valid]
       # print(f"[DEBUG] flat_gc_distances ({len(flat_gc_distances)}): "
         #     f"{flat_gc_distances[:5]}{'...' if len(flat_gc_distances)>5 else ''}")
        #print(f"[DEBUG] flat_gc_lengths    ({len(flat_gc_lengths)}):  "
          #    f"{flat_gc_lengths[:5]}{'...' if len(flat_gc_lengths)>5 else ''}")
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


        # Combine B's calculated from distant genes, and genes within the region!
    safe_flank_B = np.concatenate((np.ones(chunk_end+1 - chunk_start, dtype=float), B)) # Add an array of flank_B where all sites are B = 1, to account for sites with no flanking genes
    new_flanking_mask = np.concatenate((np.ones((1, chunk_end+1 - chunk_start), dtype=bool), mask), axis=0)
    unique_indices, inverse_indices = np.unique(np.where(new_flanking_mask)[1], return_inverse=True)
    aggregated_B = np.ones_like(np.ones_like(np.arange(chunk_start,chunk_end+1), dtype=np.float64), dtype=np.float64)
    np.multiply.at(aggregated_B, inverse_indices, safe_flank_B) # Multiplicative sum of B calculated at a given site from multiple elements

    if not quiet: print(f"Chunk {chunk_idx}: Final B_noninterfering = {aggregated_B}")
    return aggregated_B

## WON'T WORK FOR UNLINKED B, OR FOR SITES considering B from an interfering region, that aren't themselves in the interfering region
