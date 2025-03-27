import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.collections import LineCollection

def plotB(b_values_input, caller, output_path, silent, gene_ranges=None):
    if not silent: 
        print('====== P L O T T I N G . . . =======================')
    
    # Set the font family with a fallback list.
    mpl.rcParams['font.family'] = ['Helvetica', 'DejaVu Sans', 'Arial']
    
    # Check if seaborn whitegrid style is available.
    if 'seaborn-v0_8-whitegrid' in plt.style.available:
        plt.style.use('seaborn-v0_8-whitegrid')
    else:
        print("'seaborn-whitegrid' style not found, using 'ggplot' as fallback.")
        plt.style.use('ggplot')
        # Override ggplot's defaults and add grids.
        mpl.rcParams['axes.facecolor'] = 'white'
        mpl.rcParams['figure.facecolor'] = 'white'
        mpl.rcParams['grid.color'] = 'grey'
        mpl.rcParams['grid.linestyle'] = '--'
        mpl.rcParams['grid.linewidth'] = 0.5

    # Ensure axes have a black edge.
    mpl.rcParams['axes.edgecolor'] = 'black'
    
    # Create the plot.
    fig, ax = plt.subplots(figsize=(10, 6))

    # Extract x and y values.
    if caller == "genome":
        x = b_values_input['Position']
        y = b_values_input['B']
    elif caller == "region":
        x = b_values_input[0]
        y = b_values_input[1]

    # Downsample the data if there are too many points.
    max_points = 10000  # adjust this threshold as needed
    if len(x) > max_points:
        indices = np.linspace(0, len(x) - 1, max_points).astype(int)
        x = x[indices]
        y = y[indices]

    # Plot the entire B line as a solid blue line.
    ax.plot(x, y, color='blue', lw=1.5, alpha=0.8)
    ax.set_xlim(x.min() - 1, x.max())
    
    ax.set_ylabel('Expected diversity relative to neutral evolution (B)', fontsize=13)
    if caller == "genome":
        ax.set_title('B recovery across chromosomal region', fontsize=15, fontweight='bold')
        ax.set_xlabel('Chromosomal position (bp)', fontsize=13)
    elif caller == "region":
        ax.set_xlabel('Distance from single selected element of size', fontsize=13)
        ax.set_title('B recovery from single element', fontsize=15, fontweight='bold')

    ax.tick_params(axis='both', which='major', labelsize=10)

    # If gene annotations are provided, add them using LineCollection for efficiency.
    if gene_ranges is not None and len(gene_ranges) > 0:
        # Get current y-limits so we can draw below the axis.
        ymin, ymax = ax.get_ylim()
        bar_y = ymin - (ymax - ymin) * 0.05  # 5% below the axis
        
        # Expand the y-limits to make space for the gene bars.
        ax.set_ylim(bar_y, ymax)
        
        # Create segments for horizontal bars.
        segments = [((start, bar_y), (end, bar_y)) for start, end in gene_ranges]
        lc = LineCollection(segments, colors='black', linewidths=30)
        ax.add_collection(lc)
        
        # When caller is "genome", overlay black line segments for gene regions.
        if caller == "genome":
            # Create a boolean mask for x-values that fall within any gene interval.
            gene_mask = np.zeros_like(x, dtype=bool)
            for start, end in gene_ranges:
                gene_mask |= ((x >= start) & (x <= end))
            idx = np.where(gene_mask)[0]
            if len(idx) > 0:
                # Split the indices into contiguous segments.
                segs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
                # Build segments (as arrays of coordinate pairs) for the overlay.
                black_segments = [np.column_stack((x[seg], y[seg])) for seg in segs if len(seg) > 1]
                if black_segments:
                    lc_black = LineCollection(black_segments, colors='black', linewidths=1.5)
                    ax.add_collection(lc_black)

    # Format the x-axis ticks.
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda x_val, pos: f"{int(x_val)} bp" if x_val < 1000 
            else (f"{x_val/1e6:.2f} Mb" if x_val >= 1000000 else f"{int(x_val/1000)} kb")
        )
    )
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")
