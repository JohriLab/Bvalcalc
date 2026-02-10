import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib import gridspec  # for rec rate strip
from matplotlib.cm import ScalarMappable

def plotB(b_values_input, caller, output_path, quiet, gene_ranges=None, neutral_only=False, rec_rates=None, chunk_size=None, no_title=False):
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
    # Set all text and axes to black
    mpl.rcParams.update({
        'axes.edgecolor': 'black',
        'text.color': 'black',
        'axes.labelcolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black'
    })

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
    ax.set_ylabel(r'$\boldsymbol{B}$', fontsize=17, rotation=0, ha='right')
    if caller == "chromosome":
        if not no_title:
            ax.set_title(chrom, fontsize=15, fontweight='bold')
        if rec_rates is not None:
            ax.set_xlabel('Chromosomal position (bp)', fontsize=14, labelpad=38)
        else:
            ax.set_xlabel('Chromosomal position (bp)', fontsize=14, labelpad=5)
    else:
        ax.set_xlabel('Distance from single selected element of size 10 kb', fontsize=14)
        # ax.set_title('B recovery from single element', fontsize=15, fontweight='bold')

    ax.tick_params(axis='both', which='major', labelsize=11)

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
    im = None
    if rec_rates is not None and caller == "chromosome":
        ax_rec = fig.add_subplot(gs[1], sharex=ax)
        ax_rec.set_yticks([])
        ax_rec.tick_params(axis='x', which='major', labelsize=11)
        rec_img = np.expand_dims(rec_rates, axis=0)
        min_pos = positions.min()
        extent = [min_pos, min_pos + len(rec_rates) * chunk_size, 0, 1]
        im = ax_rec.imshow(rec_img, aspect='auto', extent=extent, cmap=magenta_map, origin='lower', zorder=2, vmin=0, vmax=np.max(rec_rates))
        ax_rec.set_frame_on(False)
        plt.setp(ax.get_xticklabels(), visible=False)

    # Final layout and save
    if rec_rates is None:
        plt.tight_layout()
    else:
        # Optimize layout and fit rec strip in gap between x-axis and xlabel
        plt.tight_layout()
        fig.canvas.draw()
        
        # Calculate bottom margin to fit rec strip and colorbar
        xlabel_bbox = ax.xaxis.label.get_window_extent().transformed(fig.transFigure.inverted())
        rec_strip_pos = gs[1].get_position(fig)
        bottom_for_strip = max(0.05, xlabel_bbox.y0 - rec_strip_pos.height - 0.005)
        bottom_for_colorbar = rec_strip_pos.height + 0.10  # Space for colorbar + labels
        bottom_value = max(bottom_for_strip, bottom_for_colorbar)
        fig.subplots_adjust(bottom=bottom_value, hspace=0.01)
        
        # Add colorbar at bottom
        if im is not None:
            fig.canvas.draw()
            ax_pos = ax.get_position()
            rec_pos = ax_rec.get_position()
            
            # Calculate max rec rate within the plotted region
            xlim = ax.get_xlim()
            # Find which rec_rates chunks fall within the visible xlim
            chunk_starts = np.arange(min_pos, min_pos + len(rec_rates) * chunk_size, chunk_size)
            chunk_ends = chunk_starts + chunk_size
            # Chunks that overlap with visible region
            visible_mask = (chunk_ends > xlim[0]) & (chunk_starts < xlim[1])
            max_rec_rate = np.max(rec_rates[visible_mask]) if visible_mask.any() else np.max(rec_rates)
            max_all_rec_rates = np.max(rec_rates)
            
            cbar_width = 0.144  # 20% wider than original 0.12
            cbar_height = rec_pos.height
            cbar_x = ax_pos.x1 - cbar_width
            cbar_y = 0.06
            
            cax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
            # Create a truncated colormap that only uses the portion from 0 to max_rec_rate
            # This fraction of the colormap will be stretched to fill the full colorbar
            colormap_fraction = max_rec_rate / max_all_rec_rates if max_all_rec_rates > 0 else 1.0
            # Create a colormap that samples only the first portion of the original
            truncated_colors = magenta_map(np.linspace(0, colormap_fraction, 256))
            truncated_cmap = mpl.colors.LinearSegmentedColormap.from_list('truncated', truncated_colors, N=256)
            # Create a new mappable with the visible range limits for the colorbar
            sm = ScalarMappable(cmap=truncated_cmap, norm=mpl.colors.Normalize(vmin=0, vmax=max_rec_rate))
            sm.set_array([])
            cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
            cbar.outline.set_visible(False)
            cbar.set_ticks([0, max_rec_rate])
            cbar.set_ticklabels(['0', f'{max_rec_rate:.2g} x $\\boldsymbol{{r}}$'])
            cbar.ax.tick_params(bottom=True, labelbottom=True, labelsize=11, 
                               color='black', labelcolor='black', length=4, width=1)
            cbar.ax.text(0.5, -0.3, 'CO rate', fontsize=12, ha='center', va='top', 
                        transform=cbar.ax.transAxes, color='black')

    plt.savefig(output_path, dpi=300)
    if caller == "chromosome" and not no_title and not quiet:
        print("To hide plot title, add --no_title.")
    print(f"Plot saved to {output_path}")
