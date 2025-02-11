import matplotlib.pyplot as plt

def plotBasic(region_output, output_file='../../bin/plot.png'):
    print('====== P L O T T I N G . . . =======================')
    
    plt.figure(figsize=(10, 6))
    plt.plot(region_output, color='blue', lw=1.5, label='B Recovery Curve')
    plt.xlabel('Distance from selected element (bp)')
    plt.ylabel('Relative diversity (B)')
    plt.title('B recovery from single element')
    plt.legend()
    plt.grid(True)
    
    # Save the plot to the file specified by output_file
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")
