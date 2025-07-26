import numpy as np
from Bvalcalc.core.calculateB import calculateB_linear, calculateB_recmap

def calc_B_in_genes(
    chunk_size, num_chunks, precise_chunks,
    precise_blockstart, precise_blockend,
    chunk_start, chunk_end,
    physical_lengths, precise_region_start, chunk_num,
    rec_rate_per_chunk=None, gc_rate_per_chunk=None,
    rec_lengths=None, gc_lengths=None
):
    genes_in_this_chunk_mask = np.logical_and(
        precise_blockstart < chunk_end,
        precise_blockend >= chunk_start
    )
    this_chunk_blockstart = precise_blockstart[genes_in_this_chunk_mask]
    this_chunk_blockend = precise_blockend[genes_in_this_chunk_mask]

    if len(this_chunk_blockstart) == 0:
        return np.ones(chunk_end - chunk_start, dtype=np.float64)

    this_chunk_blockstart_inchunk = np.clip(this_chunk_blockstart, chunk_start, chunk_end - 1)
    this_chunk_blockend_inchunk = np.clip(this_chunk_blockend, chunk_start, chunk_end - 1)

    within_gene_B = np.ones(chunk_end - chunk_start, dtype=np.float64) # Initialize array
    chunk_starts = precise_region_start + np.arange(num_chunks + 1) * chunk_size
    this_chunk_idx = np.searchsorted(chunk_starts, chunk_start)

    # Pre-slice once
    masked_phys_lengths = physical_lengths[genes_in_this_chunk_mask]
    masked_rec_lengths = rec_lengths[genes_in_this_chunk_mask] if rec_lengths is not None else None
    masked_gc_lengths = gc_lengths[genes_in_this_chunk_mask] if gc_lengths is not None else None

    # Precompute rate slices
    precise_slice_start = max(0, chunk_num - precise_chunks)
    precise_slice_end = min(num_chunks, chunk_num + precise_chunks + 1)

    if rec_rate_per_chunk is not None:
        precise_rates = rec_rate_per_chunk[precise_slice_start:precise_slice_end]
        rec_bp_to_element = precise_rates[this_chunk_idx]

    if gc_rate_per_chunk is not None:
        precise_gc_rates = gc_rate_per_chunk[precise_slice_start:precise_slice_end]
        gc_bp_to_element = precise_gc_rates[this_chunk_idx]

    for gene_idx in range(len(this_chunk_blockstart_inchunk)): # Ragged array needs for loop
        start = this_chunk_blockstart_inchunk[gene_idx]
        end = this_chunk_blockend_inchunk[gene_idx] + 1
        gpos_in_chunk = np.arange(start, end)
        gene_sites = gpos_in_chunk - chunk_start

        left_block_lengths = gpos_in_chunk - this_chunk_blockstart[gene_idx]
        right_block_lengths = this_chunk_blockend[gene_idx] - gpos_in_chunk
        focal_block_physlength = masked_phys_lengths[gene_idx]

        if rec_rate_per_chunk is not None and gc_rate_per_chunk is None:
            focal_block_reclength = masked_rec_lengths[gene_idx]
            scale = focal_block_reclength / focal_block_physlength
            left_chunk_reclengths = left_block_lengths * scale
            right_chunk_reclengths = right_block_lengths * scale

            left_block_B = calculateB_recmap(1, left_block_lengths, rec_bp_to_element, left_chunk_reclengths)
            right_block_B = calculateB_recmap(1, right_block_lengths, rec_bp_to_element, right_chunk_reclengths)

        elif rec_rate_per_chunk is not None and gc_rate_per_chunk is not None:
            focal_block_reclength = masked_rec_lengths[gene_idx]
            focal_block_gclength = masked_gc_lengths[gene_idx]
            # Fix needed here for RuntimeWarning: invalid value encountered in scalar divide
            rec_scale = focal_block_reclength / focal_block_physlength
            gc_scale = focal_block_gclength / focal_block_physlength

            left_chunk_reclengths = left_block_lengths * rec_scale
            right_chunk_reclengths = right_block_lengths * rec_scale

            left_chunk_gclengths = left_block_lengths * gc_scale
            right_chunk_gclengths = right_block_lengths * gc_scale

            left_block_B = calculateB_recmap(
                1, left_block_lengths, rec_bp_to_element, left_chunk_reclengths,
                gc_bp_to_element, left_chunk_gclengths
            )
            right_block_B = calculateB_recmap(
                1, right_block_lengths, rec_bp_to_element, right_chunk_reclengths,
                gc_bp_to_element, right_chunk_gclengths
            )

        elif rec_rate_per_chunk is None and gc_rate_per_chunk is not None:
            focal_block_gclength = masked_gc_lengths[gene_idx]
            gc_scale = focal_block_gclength / focal_block_physlength

            left_chunk_gclengths = left_block_lengths * gc_scale
            right_chunk_gclengths = right_block_lengths * gc_scale

            left_block_B = calculateB_recmap(1, left_block_lengths, None, None, gc_bp_to_element, left_chunk_gclengths)
            right_block_B = calculateB_recmap(1, right_block_lengths, None, None, gc_bp_to_element, right_chunk_gclengths)

        else:
            left_block_B = calculateB_linear(1, left_block_lengths)
            right_block_B = calculateB_linear(1, right_block_lengths)

        np.multiply.at(within_gene_B, gene_sites, left_block_B)
        np.multiply.at(within_gene_B, gene_sites, right_block_B)

    return within_gene_B
