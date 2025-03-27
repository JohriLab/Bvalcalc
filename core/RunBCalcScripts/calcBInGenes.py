import numpy as np
from core.calculateB import calculateB_linear, calculateB_recmap

def calcBInGenes(chunk_size, num_chunks, precise_chunks, precise_blockstart, precise_blockend, chunk_start, chunk_end, physical_lengths, precise_region_start, chunk_num, rec_rate_per_chunk = None, gc_rate_per_chunk = None, rec_lengths = None, gc_lengths = None):
    
    genes_in_this_chunk_mask = np.logical_and(precise_blockstart < chunk_end, precise_blockend >= chunk_start)
    this_chunk_blockstart = precise_blockstart[genes_in_this_chunk_mask]
    this_chunk_blockend = precise_blockend[genes_in_this_chunk_mask]

    if len(this_chunk_blockstart) == 0:
        return np.ones_like(np.arange(chunk_start, chunk_end), dtype=np.float64)

    this_chunk_blockstart_inchunk = np.clip(this_chunk_blockstart, a_min=chunk_start, a_max=chunk_end-1)
    this_chunk_blockend_inchunk = np.clip(this_chunk_blockend, a_min=chunk_start, a_max=chunk_end-1)

    agg_gene_B = np.ones_like(np.arange(chunk_start,chunk_end), dtype=np.float64)
    chunk_starts = precise_region_start + np.arange(0, num_chunks + 1) * chunk_size
    this_chunk_idx = np.searchsorted(chunk_starts, chunk_start)  # faster than np.where

    # Pre-slice these once
    masked_rec_lengths = rec_lengths[genes_in_this_chunk_mask] if rec_lengths is not None else None
    masked_gc_lengths = gc_lengths[genes_in_this_chunk_mask] if gc_lengths is not None else None
    masked_phys_lengths = physical_lengths[genes_in_this_chunk_mask]

    for gene_idx in np.arange(len(this_chunk_blockstart_inchunk)):
        gene_blockstart = this_chunk_blockstart[gene_idx]
        gene_blockend = this_chunk_blockend[gene_idx]
        gpos_in_chunk = np.arange(this_chunk_blockstart_inchunk[gene_idx], this_chunk_blockend_inchunk[gene_idx]+1)
        left_block_lengths =  gpos_in_chunk - gene_blockstart
        right_block_lengths = gene_blockend - gpos_in_chunk

        if rec_rate_per_chunk is not None and gc_rate_per_chunk is None:
            focal_block_reclength = masked_rec_lengths[gene_idx]
            focal_block_physlength = masked_phys_lengths[gene_idx]
            left_chunk_reclengths = left_block_lengths * (focal_block_reclength / focal_block_physlength)
            right_chunk_reclengths = right_block_lengths * (focal_block_reclength / focal_block_physlength)

            precise_rates = rec_rate_per_chunk[np.maximum(0, chunk_num - precise_chunks):np.minimum(num_chunks, chunk_num + precise_chunks + 1)]
            rec_bp_to_element = precise_rates[this_chunk_idx]

            left_block_B = calculateB_recmap(1, left_block_lengths, rec_bp_to_element, left_chunk_reclengths)
            right_block_B = calculateB_recmap(1, right_block_lengths, rec_bp_to_element, right_chunk_reclengths)

        elif rec_rate_per_chunk is not None and gc_rate_per_chunk is not None:
            focal_block_reclength = masked_rec_lengths[gene_idx]
            focal_block_physlength = masked_phys_lengths[gene_idx]
            left_chunk_reclengths = left_block_lengths * (focal_block_reclength / focal_block_physlength)
            right_chunk_reclengths = right_block_lengths * (focal_block_reclength / focal_block_physlength)

            precise_rates = rec_rate_per_chunk[np.maximum(0, chunk_num - precise_chunks):np.minimum(num_chunks, chunk_num + precise_chunks + 1)]
            rec_bp_to_element = precise_rates[this_chunk_idx]

            focal_block_gclength = masked_gc_lengths[gene_idx]
            left_chunk_gclengths = left_block_lengths * (focal_block_gclength / focal_block_physlength)
            right_chunk_gclengths = right_block_lengths * (focal_block_gclength / focal_block_physlength)

            precise_gc_rates = gc_rate_per_chunk[np.maximum(0, chunk_num - precise_chunks):np.minimum(num_chunks, chunk_num + precise_chunks + 1)]
            gc_bp_to_element = precise_gc_rates[this_chunk_idx]

            left_block_B = calculateB_recmap(1, left_block_lengths, rec_bp_to_element, left_chunk_reclengths, gc_bp_to_element, left_chunk_gclengths)
            right_block_B = calculateB_recmap(1, right_block_lengths, rec_bp_to_element, right_chunk_reclengths, gc_bp_to_element, right_chunk_gclengths)

        elif rec_rate_per_chunk is None and gc_rate_per_chunk is not None:
            focal_block_gclength = masked_gc_lengths[gene_idx]
            focal_block_physlength = masked_phys_lengths[gene_idx]
            left_chunk_gclengths = left_block_lengths * (focal_block_gclength / focal_block_physlength)
            right_chunk_gclengths = right_block_lengths * (focal_block_gclength / focal_block_physlength)

            precise_gc_rates = gc_rate_per_chunk[np.maximum(0, chunk_num - precise_chunks):np.minimum(num_chunks, chunk_num + precise_chunks + 1)]
            gc_bp_to_element = precise_gc_rates[this_chunk_idx]

            left_block_B = calculateB_recmap(1, left_block_lengths, None, None, gc_bp_to_element, left_chunk_gclengths)
            right_block_B = calculateB_recmap(1, right_block_lengths, None, None, gc_bp_to_element, right_chunk_gclengths)

        else:
            left_block_B = calculateB_linear(1, left_block_lengths)
            right_block_B = calculateB_linear(1, right_block_lengths)

        gene_sites = gpos_in_chunk - chunk_start
        np.multiply.at(agg_gene_B, gene_sites, left_block_B)
        np.multiply.at(agg_gene_B, gene_sites, right_block_B)

    return agg_gene_B
