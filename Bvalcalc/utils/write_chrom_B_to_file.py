import numpy as np
import sys
from .header_utils import generate_header, HeaderInfo, write_headers_to_file

def write_chrom_B_to_file(out,
                          output_data,
                          quiet,
                          hri_starts=None,
                          hri_ends=None,
                          binsize=None,
                          calc_end=None,
                          header_info=None,
                          write_header=True,
                          no_header=False):
    """
    Write Chromosome,Start,B to CSV.
    If hri_starts/hri_ends provided, append ' to B for rows whose bin overlaps any HRI region.
    - output_data: structured array with fields ['Chromosome','Start','B']
    - hri_starts, hri_ends: 1D int arrays of equal length (inclusive coordinates)
    - binsize: int, bin length (last bin truncated)
    - calc_end: int, max coordinate to cap last-bin end (optional but recommended)
    - header_info: HeaderInfo object with header information (optional)
    - write_header: Whether to write headers (default True)
    - no_header: If True, skip writing comment headers (version, command, etc.)
    """

    # Write headers if requested and not disabled
    if write_header and not no_header and header_info:
        header_lines = generate_header(header_info)
        write_headers_to_file(out, header_lines, 'w')
    elif write_header and not no_header and not header_info:
        # Write basic header if no header_info provided
        try:
            from Bvalcalc import __version__
            version_str = __version__
        except ImportError:
            version_str = "(version not found)"
        with open(out, 'w') as f:
            f.write(f"# Bvalcalc v{version_str}\n")
            f.write("# Format: Chromosome,Start,B\n")

    # Fast path: no HRI spans or no binsize provided -> write floats directly
    need_mark = (
        hri_starts is not None and hri_ends is not None and
        len(hri_starts) > 0 and len(hri_ends) > 0 and
        binsize is not None
    )

    if not need_mark:
        # Always use append mode (header either written here or by genomeBcalc)
        mode = 'a'
        with open(out, mode) as f:
            np.savetxt(f, output_data, delimiter=",", fmt="%s,%d,%.6f", comments="")
        return

    # Build bin start/end coordinates
    bin_starts = output_data['Start'].astype(np.int64)
    if calc_end is None:
        # Fallback: cap ends at start + binsize - 1
        bin_ends = bin_starts + (int(binsize) - 1)
    else:
        bin_ends = np.minimum(bin_starts + (int(binsize) - 1), int(calc_end)).astype(np.int64)

    # Sort HRI intervals for binary searches
    s = np.sort(np.asarray(hri_starts, dtype=np.int64))
    e = np.sort(np.asarray(hri_ends,   dtype=np.int64))

    # Vectorized interval overlap test:
    # overlap if count(starts <= bin_end) - count(ends < bin_start) > 0
    starts_le_end = np.searchsorted(s, bin_ends, side='right')
    ends_lt_start = np.searchsorted(e, bin_starts - 1, side='right')
    overlaps = (starts_le_end - ends_lt_start) > 0

    # Prepare per-row B strings with optional '*'
    B = output_data['B'].astype(float)
    B_as_str = np.array([f"{v:.6f}" for v in B], dtype="<U16")
    if overlaps.any():
        B_as_str[overlaps] = np.char.add(B_as_str[overlaps], "'")

    # Compose rows (don’t mutate output_data dtype)
    rows = np.empty((output_data.shape[0], 3), dtype=object)
    rows[:, 0] = output_data['Chromosome']
    rows[:, 1] = output_data['Start'].astype(int)
    rows[:, 2] = B_as_str

    if not quiet and hri_starts is not None: 
        print(f"B values calculated for HRI regions (B') is indicated with ' at end of line")

    # Always use append mode (header either written here or by genomeBcalc)
    mode = 'a'
    with open(out, mode) as f:
        np.savetxt(f, rows, delimiter=",", fmt="%s,%d,%s", comments="")

