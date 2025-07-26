import numpy as np
import sys

def bin_outputs(b_values: np.ndarray,
                positions: np.ndarray,
                binsize: int
               ) -> tuple[np.ndarray, np.ndarray]:
    """
    Collapse b_values & positions into consecutive bins of length `binsize`.
    Returns
    -------
    binned_b_values : float array, shape=(ceil(N/binsize),)
        Mean of b_values in each bin.
    binned_positions : int array, same length as binned_b_values
        The first position in each bin (i.e. positions[0], positions[binsize], …).
    """
    if b_values.shape[0] != positions.shape[0]:
        print("Shapes", b_values.shape[0], positions.shape[0])
        raise ValueError("`b_values` and `positions` must be the same length")

    N = b_values.shape[0]
    # bin start indices: 0, binsize, 2*binsize, …
    bin_idx = np.arange(0, N, binsize, dtype=int)

    # fast C‐loop to sum each bin
    sums = np.add.reduceat(b_values, bin_idx)

    # how many items in each bin (last bin may be shorter)
    counts = np.diff(np.append(bin_idx, N))

    # compute mean
    binned_b_values = sums / counts

    # record the first genomic position in each bin
    binned_positions = positions[bin_idx]

    return binned_b_values, binned_positions
