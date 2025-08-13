import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib import gridspec  # for rec rate strip

def plotB(b_values_input, caller, output_path, quiet, gene_ranges=None, neutral_only=False, rec_rates=None, chunk_size=None):
    if not quiet:
        print('====== P L O T T I N G . . . =======================')

    # Configure fonts and styles
    mpl.rcParams['font.family'] = ['Helvetica', 'DejaVu Sans', 'Arial']
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

    # Create figure and axis based on rec_rates
    if rec_rates is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        height_ratios = [10, 0.5]
        fig = plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=height_ratios)
        ax = fig.add_subplot(gs[0])

    # Main plotting logic
    if caller == "chromosome":
        positions = b_values_input['Start']
        b_vals = b_values_input['B']
        chrom = b_values_input['Chromosome'][0] if 'Chromosome' in b_values_input.dtype.names else 'unknown'

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
            for xs, ys in zip(x_segments, y_segments):
                if len(xs) > 1:
                    ax.plot(xs, ys, color='blue', lw=1.5, alpha=0.8)

            if len(x) > 0:
                ax.set_xlim(x.min() - 1, x.max())
        else:
            x = positions
            y = b_vals
            max_points = 10000
            if len(x) > max_points:
                idx = np.linspace(0, len(x) - 1, max_points).astype(int)
                x = x[idx]
                y = y[idx]
            ax.plot(x, y, color='blue', lw=1.5, alpha=0.8)
            ax.set_xlim(x.min() - 1, x.max())

    elif caller == "gene":
        x = b_values_input[:, 0]
        y = b_values_input[:, 1]
        ax.plot(x, y, color='blue', lw=1.5, alpha=0.8)
        ax.set_xlim(x.min() - 1, x.max())

    # Labels and title
    ax.set_ylabel('Expected diversity relative to neutral evolution (B)', fontsize=13)
    if caller == "chromosome":
        ax.set_title(f'B for chromosome {chrom} ({positions.min()}–{positions.max()} bp)', fontsize=15, fontweight='bold')
        if rec_rates is not None:
            ax.set_xlabel('Chromosomal position (bp)', fontsize=13, labelpad=40)
        else:
            ax.set_xlabel('Chromosomal position (bp)', fontsize=13)
    else:
        ax.set_xlabel('Distance from single selected element of size', fontsize=13)
        ax.set_title('B recovery from single element', fontsize=15, fontweight='bold')

    ax.tick_params(axis='both', which='major', labelsize=10)

    # Gene-range bars
    if gene_ranges is not None and len(gene_ranges) > 0:
        if caller == "chromosome":
            ymin, ymax = ax.get_ylim()
            bar_y = ymin - (ymax - ymin) * 0.05

            # compute the current x-range
            xlim = ax.get_xlim()
            xrange = xlim[1] - xlim[0]
            fig_width_in_pixels = fig.bbox.width
            min_width = xrange / fig_width_in_pixels

            # convert start/end to int for plotting, ensuring min width
            segments = []
            for _, start, end in gene_ranges:
                s = float(start)
                e = float(end)
                if (e - s) < min_width:
                    e = s + min_width
                segments.append(((s, bar_y), (e, bar_y)))

            ax.set_ylim(bar_y, ymax)
            ax.add_collection(LineCollection(segments, colors='black', linewidths=30))

            if not neutral_only:
                gene_mask = np.zeros_like(x, dtype=bool)
                for _, start, end in gene_ranges:
                    start = int(start)
                    end = int(end)
                    gene_mask |= (x >= start) & (x <= end)
                idx = np.where(gene_mask)[0]
                splits = np.where(np.diff(idx) != 1)[0] + 1
                for seg in np.split(idx, splits):
                    if len(seg) > 1:
                        coords = np.column_stack((x[seg], y[seg]))
                        ax.add_collection(LineCollection([coords], colors='black', linewidths=1.5))
        else:
            # For caller != chromosome, keep original behavior (no min width enforcement)
            ymin, ymax = ax.get_ylim()
            bar_y = ymin - (ymax - ymin) * 0.05
            ax.set_ylim(bar_y, ymax)
            segments = [((int(start), bar_y), (int(end), bar_y)) for _, start, end in gene_ranges]
            ax.add_collection(LineCollection(segments, colors='black', linewidths=30))

    # X-axis formatting
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda value, pos:
            f"{int(value)} bp" if value < 1e3 else (
                f"{value/1e6:.2f} Mb" if value >= 1e6 else f"{int(value/1e3)} kb"
            )
        )
    )

    from matplotlib.colors import LinearSegmentedColormap
    magenta_map = LinearSegmentedColormap.from_list("custom_magenta", ["white", "#C54B8C"])
    
    # Recombination rate strip
    if rec_rates is not None and caller == "chromosome":
        ax_rec = fig.add_subplot(gs[1], sharex=ax)
        ax_rec.set_yticks([])
        ax_rec.tick_params(axis='x', which='major', labelsize=9)
        rec_img = np.expand_dims(rec_rates, axis=0)
        min_pos = positions.min()
        extent = [min_pos, min_pos + len(rec_rates) * chunk_size, 0, 1]
        ax_rec.imshow(rec_img, aspect='auto', extent=extent, cmap=magenta_map, origin='lower', zorder=2, vmin=0, vmax=np.max(rec_rates))
        ax_rec.set_frame_on(False)
        plt.setp(ax.get_xticklabels(), visible=False)

    # Final layout and save
    if rec_rates is None:
        try:
            plt.savefig(output_path, dpi=300)
        except Exception as e:           # you could catch TimeoutError specifically if you prefer
            import sys
            print(f"Warning: could not save plot to {output_path!r}: {e}", file=sys.stderr)
    else:
        fig.subplots_adjust(hspace=0.01, bottom=0.12)

    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")
