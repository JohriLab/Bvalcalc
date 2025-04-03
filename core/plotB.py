import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.collections import LineCollection

def plotB(b_values_input, caller, output_path, quiet, gene_ranges=None, neutral_only=False):
    if not quiet: 
        print('====== P L O T T I N G . . . =======================')

    # Set the font family with a fallback list.
    mpl.rcParams['font.family'] = ['Helvetica', 'DejaVu Sans', 'Arial']

    # Style
    if 'seaborn-v0_8-whitegrid' in plt.style.available:
        plt.style.use('seaborn-v0_8-whitegrid')
    else:
        print("'seaborn-whitegrid' style not found, using 'ggplot' as fallback.")
        plt.style.use('ggplot')
        mpl.rcParams['axes.facecolor'] = 'white'
        mpl.rcParams['figure.facecolor'] = 'white'
        mpl.rcParams['grid.color'] = 'grey'
        mpl.rcParams['grid.linestyle'] = '--'
        mpl.rcParams['grid.linewidth'] = 0.5

    mpl.rcParams['axes.edgecolor'] = 'black'

    fig, ax = plt.subplots(figsize=(10, 6))

    # Genome mode: extract and optionally filter/segment
    if caller == "genome":
        positions = b_values_input['Position']
        b_vals = b_values_input['B']

        if neutral_only:
            conserved = b_values_input['Conserved']
            is_bytes = conserved.dtype.kind == 'S'
            neutral_mask = conserved == b'N' if is_bytes else conserved == 'N'

            x = positions[neutral_mask]
            y = b_vals[neutral_mask]

            sort_idx = np.argsort(x)
            x = x[sort_idx]
            y = y[sort_idx]

            diffs = np.diff(x)
            split_indices = np.where(diffs != 1)[0] + 1
            x_segments = np.split(x, split_indices)
            y_segments = np.split(y, split_indices)

            print(f"Plotting {len(x)} neutral positions in {len(x_segments)} segments.")
            for i, (x_seg, y_seg) in enumerate(zip(x_segments, y_segments)):
                if len(x_seg) > 1:
                    ax.plot(x_seg, y_seg, color='blue', lw=1.5, alpha=0.8)

            if len(x) > 0:
                ax.set_xlim(x.min() - 1, x.max())

        else:
            x = positions
            y = b_vals

            max_points = 10000
            if len(x) > max_points:
                indices = np.linspace(0, len(x) - 1, max_points).astype(int)
                x = x[indices]
                y = y[indices]

            ax.plot(x, y, color='blue', lw=1.5, alpha=0.8)
            ax.set_xlim(x.min() - 1, x.max())

    # Region mode
    elif caller == "region":
        x = b_values_input[:, 0]  # column 0 = positions
        y = b_values_input[:, 1]  # column 1 = B values

        ax.plot(x, y, color='blue', lw=1.5, alpha=0.8)
        ax.set_xlim(x.min() - 1, x.max())

    # Labels and title
    ax.set_ylabel('Expected diversity relative to neutral evolution (B)', fontsize=13)

    if caller == "genome":
        ax.set_title('B recovery across chromosomal region', fontsize=15, fontweight='bold')
        ax.set_xlabel('Chromosomal position (bp)', fontsize=13)
    elif caller == "region":
        ax.set_xlabel('Distance from single selected element of size', fontsize=13)
        ax.set_title('B recovery from single element', fontsize=15, fontweight='bold')

    ax.tick_params(axis='both', which='major', labelsize=10)

    # Gene bars
    if gene_ranges is not None and len(gene_ranges) > 0:
        ymin, ymax = ax.get_ylim()
        bar_y = ymin - (ymax - ymin) * 0.05
        ax.set_ylim(bar_y, ymax)

        segments = [((start, bar_y), (end, bar_y)) for start, end in gene_ranges]
        lc = LineCollection(segments, colors='black', linewidths=30)
        ax.add_collection(lc)

        if caller == "genome" and not neutral_only:
            gene_mask = np.zeros_like(x, dtype=bool)
            for start, end in gene_ranges:
                gene_mask |= ((x >= start) & (x <= end))
            idx = np.where(gene_mask)[0]
            if len(idx) > 0:
                segs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
                black_segments = [np.column_stack((x[seg], y[seg])) for seg in segs if len(seg) > 1]
                if black_segments:
                    lc_black = LineCollection(black_segments, colors='black', linewidths=1.5)
                    ax.add_collection(lc_black)

    # Format the x-axis ticks
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda x_val, pos: f"{int(x_val)} bp" if x_val < 1000 
            else (f"{x_val/1e6:.2f} Mb" if x_val >= 1e6 else f"{int(x_val/1e3)} kb")
        )
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")
