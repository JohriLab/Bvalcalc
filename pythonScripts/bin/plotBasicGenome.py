import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker

def plotBasicGenome(genome_output, output_file='../../bin/plotGenome.png'):
    print('====== P L O T T I N G . . . =======================')
    
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
    ax.plot(genome_output, color='black', lw=1.5)
    ax.set_xlabel('Distance from selected element', fontsize=13)
    ax.set_ylabel('Relative diversity (B)', fontsize=13)
    ax.set_title('B across chromosomal region', fontsize=15, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Format the x-axis ticks:
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda x, pos: f"{int(x)} bp" if x < 1000 else (f"{int(x/1e6)} Mb" if x >= 1000000 else f"{int(x/1000)} kb")
        )
    )
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")
