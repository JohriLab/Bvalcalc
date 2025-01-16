import numpy as np

def calculateLPerChunk(
    chunk_size,
    blockstart,
    blockend,
    chr_start,
    chr_end
):
    """
    Loop-based calculation of how much each set of blocks (flank blocks + main blocks)
    overlaps with each chunk in [chr_start, chr_end).

    Parameters
    ----------
    chunk_size : int
        Size of each chunk.
    flank_blockstart : 1D array-like
        Start positions for the flank blocks.
    flank_blockend : 1D array-like
        End positions for the flank blocks.
    blockstart : 1D array-like
        Start positions for the main blocks.
    blockend : 1D array-like
        End positions for the main blocks.
    chr_start : int
        Start of the chromosome/region of interest.
    chr_end : int
        End of the chromosome/region of interest (non-inclusive boundary).

    Returns
    -------
    chunk_sums : np.ndarray
        1D array with length = number of chunks.
        Each entry is the total overlap of all blocks with that chunk.
    """

    # Number of chunks, ensuring we cover the entire [chr_start, chr_end)
    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size  # ceiling division
    
    # Precompute the start and end of each chunk
    chunk_indices = np.arange(num_chunks)
    chunk_starts = chr_start + chunk_indices * chunk_size
    chunk_ends   = chr_start + (chunk_indices + 1) * chunk_size

    # Allocate array to accumulate overlap length per chunk
    chunk_sums = np.zeros(num_chunks, dtype=float)

    def accumulate_overlaps(starts, ends):
        """
        Accumulate overlaps for a set of blocks with the pre-defined chunks.
        """
        for s, e in zip(starts, ends):
            # Ignore blocks with no positive length or completely out of range
            if e <= s:
                continue
            if e <= chr_start or s >= chr_end:
                continue

            # Clamp block coordinates to [chr_start, chr_end)
            s_clamped = max(s, chr_start)
            e_clamped = min(e, chr_end)

            if e_clamped <= s_clamped:
                continue  # No valid overlap

            # Determine which chunks this block might span
            # (integer division for chunk indices)
            start_chunk_idx = (s_clamped - chr_start) // chunk_size
            end_chunk_idx   = (e_clamped - 1 - chr_start) // chunk_size

            # Ensure indices are in [0, num_chunks-1]
            start_chunk_idx = max(0, start_chunk_idx)
            end_chunk_idx   = min(num_chunks - 1, end_chunk_idx)

            # Walk through the relevant chunks and accumulate overlap
            for c in range(start_chunk_idx, end_chunk_idx + 1):
                overlap_start = max(s_clamped, chunk_starts[c])
                overlap_end   = min(e_clamped, chunk_ends[c])
                overlap_len   = overlap_end - overlap_start
                if overlap_len > 0:
                    chunk_sums[c] += overlap_len

    # 2) Accumulate overlap within given chunk for genes
    accumulate_overlaps(blockstart, blockend)

    return chunk_sums