import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib import gridspec  # for rec rate strip

def plotB_figures(b_values_input, caller, output_path, quiet, gene_ranges=None, neutral_only=False, rec_rates=None):
    if not quiet:
        print('====== P L O T T I N G . . . =======================ss')

    B_uncorrected = None
    B_observed = None
    legend_name_blue = "Calculated B"
    legend_name_orange = "Uncorrected B"
    legend_name_dot = "Observed B"

    # nogcBasicParams AKA Normal
    # poetry run Bvalcalc --gene --params ./tests/testparams/nogcBasicParams.py --plot /Users/jmarsh96/Desktop/Bcalc/Figures/nogcBasicParams.png
    B_observed = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/all.pi"
    legend_name_blue = "Calculated"
    legend_name_dot = "Observed (simulations)"
    title_name = 'B recovery from single element with DFE'

    # gcBasicParams AKA with GC
    # poetry run Bvalcalc --gene --params ./tests/testparams/gcBasicParams.py --plot /Users/jmarsh96/Desktop/Bcalc/Figures/gcBasicParams.png
    # B_uncorrected = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/nogcBasicParams.B"
    # B_observed = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/40kb_gc_all.pi"
    # legend_name_blue = "Calculated (with GC)"
    # legend_name_orange = "Calculated (no GC)"
    # legend_name_dot = "Observed (simulations)"
    # title_name = 'B recovery from single element with gene conversion'

    # expand_5N_1T AKA Demography
    # poetry run Bvalcalc --gene --params tests/testparams/ExpandParams_5N_1T.py --pop_change --plot /Users/jmarsh96/Desktop/Bcalc/Figures/expand_5N_1T.png
    # B_uncorrected = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/c40kb_expand_5N_1T.bvals"
    # B_observed = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/OBS_Expand_5N_1T.csv"
    # legend_name_blue = "Calculated (with demography)"
    # legend_name_orange = "Calculated (no demography)"
    # legend_name_dot = "Observed (simulations)"
    # title_name = 'B recovery from single element (5X Expansion 1N_anc generations ago)'

    # SelfParams_0.9S_0.5h AKA Selfing
    # poetry run Bvalcalc --gene --params tests/testparams/SelfParams_0.9S_0.5h.py --plot /Users/jmarsh96/Desktop/Bcalc/Figures/SelfParams_0.9S_0.5h.png
    # B_uncorrected = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/nogcBasicParams.B"
    # B_observed = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/40kb_f0.9_h0.5_all.pi"
    # legend_name_blue = "Calculated (with selfing)"
    # legend_name_orange = "Calculated (no selfing)"
    # legend_name_dot = "Observed (simulations)"
    # title_name = 'B recovery from single element with selfing (S = 0.9)'

    # Rename to reflect the actual parameter change

    if B_uncorrected is not None: load_B_uncorrected(B_uncorrected)
    if B_observed is not None: load_B_observed(B_observed)

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
        positions = b_values_input['Position']
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
                    ax.plot(xs, ys, color='blue', lw=2, zorder=1, alpha=0.8, label=legend_name_blue)#"Calculated B")

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
            ax.plot(x, y, color='blue', lw=2, zorder=1, alpha=0.8, label=legend_name_blue)#"Calculated B")
            ax.set_xlim(x.min() - 1, x.max())

    elif caller == "gene":
        x = b_values_input[:, 0]
        y = b_values_input[:, 1]
        ax.plot(x, y, color='blue', lw=2, zorder=1, alpha=0.8, label=legend_name_blue)#"Calculated B")
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
        ax.set_xlabel('Distance from single selected element of size 10 kb', fontsize=13)
        ax.set_title(title_name, fontsize=15, fontweight='bold')

        # Plot B_uncorrected if provided
    if B_uncorrected is not None:
        uncorrected_data = load_B_uncorrected(B_uncorrected)
        if uncorrected_data is not None and len(uncorrected_data) > 0:
            ax.plot(
                uncorrected_data["Distance"],
                uncorrected_data["B"],
                color='black',
                lw=2,
                alpha=0.5,
                label=legend_name_orange
            )
        # Plot B_uncorrected if provided
    if B_observed is not None:
        observed_data = load_B_observed(B_observed)
        if observed_data is not None and len(observed_data) > 0:
            ax.scatter(
                observed_data["Distance"],
                observed_data["B"],
                color='#f57616',
                s=12,
                alpha=0.9,
                label=legend_name_dot#"Observed B (simulations)"
            )


    ax.tick_params(axis='both', which='major', labelsize=10)

    # Gene-range bars
    if gene_ranges is not None and len(gene_ranges) > 0:
        ymin, ymax = ax.get_ylim()
        bar_y = ymin - (ymax - ymin) * 0.05
        ax.set_ylim(bar_y, ymax)

        # convert start/end to int for plotting
        segments = [((int(start), bar_y), (int(end), bar_y)) for _, start, end in gene_ranges]
        ax.add_collection(LineCollection(segments, colors='black', linewidths=30))

        if caller == "chromosome" and not neutral_only:
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

    # X-axis formatting
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda value, pos:
            f"{int(value)} bp" if value < 1e3 else (
                f"{value/1e6:.2f} Mb" if value >= 1e6 else f"{int(value/1e3)} kb"
            )
        )
    )

    # Recombination rate strip
    if rec_rates is not None and caller == "chromosome":
        ax_rec = fig.add_subplot(gs[1], sharex=ax)
        ax_rec.set_yticks([])
        ax_rec.tick_params(axis='x', which='major', labelsize=9)
        rec_img = np.expand_dims(rec_rates, axis=0)
        min_pos = positions.min()
        extent = [min_pos, min_pos + len(rec_rates) * 20000, 0, 1]
        ax_rec.imshow(rec_img, aspect='auto', extent=extent, cmap='Purples', origin='lower')
        ax_rec.set_frame_on(False)
        plt.setp(ax.get_xticklabels(), visible=False)

    # Final layout and save
    if rec_rates is None:
        plt.tight_layout()
    else:
        fig.subplots_adjust(hspace=0.01, bottom=0.12)
    ax.legend(loc="lower right", fontsize=10, frameon=True)


    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")

import numpy as np

def load_B_uncorrected(file_path):
    """
    Loads B_uncorrected data from a CSV file with header 'Distance,B'
    
    Parameters
    ----------
    file_path : str
        Path to the uncorrected B-values file

    Returns
    -------
    np.ndarray
        Structured array with fields 'Distance' and 'B'
    """
    try:
        data = np.genfromtxt(file_path, delimiter=",", names=True, dtype=[("Distance", float), ("B", float)])
        return data
    except Exception as e:
        print(f"Error reading B_observed file: {file_path}")
        raise e

def load_B_observed(file_path):
    """
    Loads B_observed data from a CSV file with header 'replicate,start,pi',
    averages across replicates by 'start', subtracts 9999 from 'start', and returns
    a structured array with fields 'Distance' and 'B'.
    
    Parameters
    ----------
    file_path : str
        Path to the observed B-values file
    
    Returns
    -------
    np.ndarray
        Structured array with fields 'Distance' and 'B'
    """
    try:
        raw = np.genfromtxt(file_path, delimiter=",", names=True)
        # Convert to Pandas for groupby convenience
        import pandas as pd
        df = pd.DataFrame(raw)
        grouped = df.groupby("start", as_index=False)["pi"].mean()
        grouped["Distance"] = grouped["start"] - 9999
        grouped["pi"] = grouped["pi"] / 0.012 #/ 0.00631579 # for selfing, / 0.00631579, else 0.012
        grouped.rename(columns={"pi": "B"}, inplace=True)
        result = grouped[["Distance", "B"]].to_records(index=False)
        return result
    except Exception as e:
        print(f"Error reading B_uncorrected file: {file_path}")
        raise e