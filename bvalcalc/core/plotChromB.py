import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def plotChromB(flat_b, flat_chrom, output_path):
    print('====== P L O T T I N G . . . =======================')

    mpl.rcParams['font.family'] = ['Helvetica', 'DejaVu Sans', 'Arial']

    if 'seaborn-v0_8-whitegrid' in plt.style.available:
        plt.style.use('seaborn-v0_8-whitegrid')
    else:
        print("'seaborn-whitegrid' style not found, using 'ggplot' as fallback.")
        plt.style.use('ggplot')
        mpl.rcParams['axes.facecolor']  = 'white'
        mpl.rcParams['figure.facecolor'] = 'white'
        mpl.rcParams['grid.color']      = 'grey'
        mpl.rcParams['grid.linestyle']  = '--'
        mpl.rcParams['grid.linewidth']  = 0.5

    mpl.rcParams['axes.edgecolor'] = 'black'

    fig, ax = plt.subplots(figsize=(10, 6))

    # Group B-values by chromosome
    chroms = np.unique(flat_chrom)
    data   = [flat_b[flat_chrom == chrom] for chrom in chroms]

    # Draw boxplot
    bp = ax.boxplot(
        data,
        labels=chroms,
        patch_artist=True,
        boxprops=dict(facecolor='lightblue', edgecolor='black'),
        medianprops=dict(color='red')
    )

    ax.set_xlabel('Chromosome', fontsize=13)
    ax.set_ylabel('Expected diversity relative to neutral evolution (B)', fontsize=13)
    ax.set_title('Distribution of B by Chromosome', fontsize=15, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")
