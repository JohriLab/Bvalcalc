import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib import gridspec  # for rec rate strip
import sys
from matplotlib.lines import Line2D

def plotB_figures_200kb(b_values_input, caller, output_path, quiet, gene_ranges=None, neutral_only=False, rec_rates=None):
    if not quiet:
        print('====== P L O T T I N G . . . =======================ss')

    B_uncorrected = None
    B_observed = None
    legend_name_blue = "Calculated B (intergenic)"
    legend_name_orange = "Uncorrected B"
    legend_name_dot = "Observed (simulations)"
    legend_name_black = "Calculated B (synonymous)"

    # chr_200kb AKA 200kb Genome
    # poetry run Bvalcalc --region chr_200kb:1-200000 --params tests/testparams/nogcBasicParams.py --bedgff tests/testfiles/200kb_slimtest.csv --plot /Users/jmarsh96/Desktop/Bcalc/Figures/chr_200kb.png
    B_observed = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/200kb_all.pi"
    title_name = 'B for 200 kb genome with 10 selected elements'

    # chr_200kb_recmap AKA 200kb Rec Map
    # poetry run Bvalcalc --region chr_200kb:1-200000 --params tests/testparams/nogcBasicParams.py --bedgff tests/testfiles/200kb_slimtest.csv --plot /Users/jmarsh96/Desktop/Bcalc/Figures/chr_200kb_recmap.png
    # B_observed = "/Users/jmarsh96/Desktop/Bcalc/Figures/data/200kb_recmap_all.pi"
    # title_name = 'B for 200 kb genome with recombination map'
    Genome = True

    if B_uncorrected is not None:
        load_B_uncorrected(B_uncorrected)
    if B_observed is not None:
        load_B_observed(B_observed)

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

    if rec_rates is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        height_ratios = [10, 0.5]
        fig = plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=height_ratios)
        ax = fig.add_subplot(gs[0])

    if caller == "chromosome":
        positions = b_values_input['Position']
        b_vals = b_values_input['B']
        chrom = b_values_input['Chromosome'][0] if 'Chromosome' in b_values_input.dtype.names else 'unknown'

        x = positions
        y = b_vals

        max_points = 10000
        if len(x) > max_points:
            idx = np.linspace(0, len(x) - 1, max_points).astype(int)
            x = x[idx]
            y = y[idx]
        ax.plot(x, y, color='black', lw=1.5, alpha=1, zorder=1, label=legend_name_black)
        ax.plot(x, y, color='blue', lw=1.5, alpha=1, zorder=2, label=legend_name_blue)
        ax.set_xlim(x.min() - 1, x.max())
        ax.set_ylim(0.4, 1.0)

    elif caller == "gene":
        x = b_values_input[:, 0]
        y = b_values_input[:, 1]
        ax.plot(x, y, color='blue', lw=1.5, alpha=0.8, label=legend_name_blue)
        ax.set_xlim(x.min() - 1, x.max())
        ax.set_ylim(0.4, 1.0)

    ax.set_ylabel('Expected diversity relative to neutral evolution (B)', fontsize=13)
    if caller == "chromosome":
        ax.set_title(f'{title_name}', fontsize=15, fontweight='bold')
        if rec_rates is not None:
            ax.set_xlabel('Chromosomal position (bp)', fontsize=13, labelpad=40)
        else:
            ax.set_xlabel('Chromosomal position (bp)', fontsize=13)
    else:
        ax.set_xlabel('Distance from single selected element of size 10 kb', fontsize=13)
        ax.set_title(title_name, fontsize=15, fontweight='bold')

    if B_uncorrected is not None:
        uncorrected_data = load_B_uncorrected(B_uncorrected)
        if uncorrected_data is not None and len(uncorrected_data) > 0:
            ax.plot(
                uncorrected_data["Distance"],
                uncorrected_data["B"],
                color='orange',
                lw=1.5,
                alpha=0.5,
                label=legend_name_orange
            )

    if Genome is True:
        if B_observed is not None:
            print(f"Add code")
            observed_data = load_B_observed(B_observed)
            observed_data = observed_data[observed_data["B"] > 0.001]

    else:
        if B_observed is not None:
            observed_data = load_B_observed(B_observed)

    if observed_data is not None and len(observed_data) > 0:
        ax.scatter(
            observed_data["Distance"] + 500,
            observed_data["B"],
            color='#f57616',
            s=20,
            alpha=0.9,
            label=legend_name_dot,
            zorder=3
        )

    ax.tick_params(axis='both', which='major', labelsize=10)

    if gene_ranges is not None and len(gene_ranges) > 0:
        ymin, ymax = ax.get_ylim()
        bar_y = ymin - (ymax - ymin) * 0.05
        ax.set_ylim(bar_y, ymax)

        segments = [((int(start), bar_y), (int(end), bar_y)) for _, start, end in gene_ranges]
        ax.add_collection(LineCollection(segments, colors='black', linewidths=30))

        x = np.asarray(x, dtype=float)
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
                ax.add_collection(LineCollection([coords], colors='black', linewidths=1.5, label=legend_name_black if seg is splits[0] else "_nolegend_", zorder=2))

    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda value, pos:
            f"{int(value)} bp" if value < 1e3 else (
                f"{value/1e6:.2f} Mb" if value >= 1e6 else f"{int(value/1e3)} kb"
            )
        )
    )
    from matplotlib.colors import LinearSegmentedColormap
    magenta_map = LinearSegmentedColormap.from_list("custom_magenta", ["white", "#C54B8C"])


    if rec_rates is not None and caller == "chromosome":
        ax_rec = fig.add_subplot(gs[1], sharex=ax)
        ax_rec.set_yticks([])
        ax_rec.tick_params(axis='x', which='major', labelsize=9)
        rec_img = np.expand_dims(rec_rates, axis=0)
        min_pos = positions.min()
        extent = [min_pos, min_pos + len(rec_rates) * 20000, 0, 1]
        ax_rec.imshow(rec_img, aspect='auto', extent=extent, cmap=magenta_map, origin='lower', zorder=2, vmin=0, vmax=np.max(rec_rates))
        ax_rec.set_frame_on(False)
        plt.setp(ax.get_xticklabels(), visible=False)

    if rec_rates is None:
        plt.tight_layout()
    else:
        fig.subplots_adjust(hspace=0.01, bottom=0.12)

    ax.legend(loc="lower right", bbox_to_anchor=(1, 0.04), fontsize=10, frameon=True)

    print("Observed dot at:", observed_data)
    print("Calculated B at 19999 and 20000:", b_values_input[19997:20101])

    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")

import numpy as np

def load_B_uncorrected(file_path):
    try:
        data = np.genfromtxt(file_path, delimiter=",", names=True, dtype=[("Distance", float), ("B", float)])
        return data
    except Exception as e:
        print(f"Error reading B_observed file: {file_path}")
        raise e

def load_B_observed(file_path):
    try:
        raw = np.genfromtxt(file_path, delimiter=",", names=True)
        import pandas as pd
        df = pd.DataFrame(raw)
        grouped = df.groupby("start", as_index=False)["pi"].mean()
        grouped["Distance"] = grouped["start"]
        grouped["pi"] = grouped["pi"] / 0.012
        grouped.rename(columns={"pi": "B"}, inplace=True)
        result = grouped[["Distance", "B"]].to_records(index=False)
        return result
    except Exception as e:
        print(f"Error reading B_uncorrected file: {file_path}")
        raise e
