import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker

def plotBasic(b_values_input, caller, output_path, silent):
    if not silent: print('====== P L O T T I N G . . . =======================')
    print("hiii")
    
    # Set the font family with a fallback list.
    mpl.rcParams['font.family'] = ['Helvetica', 'DejaVu Sans', 'Arial']
    
    # Check if seaborn whitegrid style is available.
    if 'seaborn-v0_8-whitegrid' in plt.style.available:
        plt.style.use('seaborn-v0_8-whitegrid')
    else:
        print("'seaborn-whitegrid' style not found, using 'ggplot' as fallback.")
        plt.style.use('ggplot')
        # Override ggplot's defaults and add grids
        mpl.rcParams['axes.facecolor'] = 'white'
        mpl.rcParams['figure.facecolor'] = 'white'
        mpl.rcParams['grid.color'] = 'grey'
        mpl.rcParams['grid.linestyle'] = '--'
        mpl.rcParams['grid.linewidth'] = 0.5

    # Ensure axes have a black edge.
    mpl.rcParams['axes.edgecolor'] = 'black'
    
    # Create the plot.
    fig, ax = plt.subplots(figsize=(10, 6))

    x = b_values_input[:, 0]  # x-values (positions)
    y = b_values_input[:, 1]  # y-values (corresponding values)

    ax.plot(x, y, color='black', lw=1.5)
    ax.set_xlim(x.min() - 1, x.max())

    ax.set_ylabel('Expected diversity relative to neutral evolution (B)', fontsize=13)
    if caller == "genome":
        ax.set_title('B recovery across chromosomal region [200 kb]', fontsize=15, fontweight='bold')
        ax.set_xlabel('Chromosomal position (bp)', fontsize=13)
    elif caller == "region":
        ax.set_xlabel('Distance from single selected element of size [40 kb]', fontsize=13)
        ax.set_title('B recovery from single element', fontsize=15, fontweight='bold')

    ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Format the x-axis ticks:
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda x, pos: f"{int(x)} bp" if x < 1000 else (f"{x/1e6:.2f} Mb" if x >= 1000000 else f"{int(x/1000)} kb")
        )
    )
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")