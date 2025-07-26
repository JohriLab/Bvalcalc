import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np

def plotChromB(flat_b, flat_chrom, output_path, quiet):
     if not quiet: print('====== P L O T T I N G . . . =======================')

     # Set overall style
     sns.set_theme(style="whitegrid", font='Helvetica')
     mpl.rcParams['axes.edgecolor'] = 'black'

     fig, ax = plt.subplots(figsize=(12, 6))

     # Unique chromosomes plus an "All" category
     chroms = list(np.unique(flat_chrom))
     chroms.append('All')

     # Extend arrays to include the combined "All" data
     all_label = np.array(['All'] * len(flat_b), dtype='<U20')
     ext_chrom = np.concatenate([flat_chrom, all_label])
     ext_b = np.concatenate([flat_b, flat_b])

     # Boxenplot including the "All" column, with hue to align palette
     # Use k_depth='tukey' to speed up computation by using Tukey's five-number summary
     sns.boxenplot(
         x=ext_chrom,
         y=ext_b,
         order=chroms,
         hue=ext_chrom,
         palette='pastel',
         showfliers=False,
         legend=False,
         k_depth='tukey',
         ax=ax
     )

     # Final plot adjustments
     ax.set_xlabel('Chromosome', fontsize=13)
     ax.set_ylabel('Expected diversity relative to neutral evolution (B)', fontsize=13)
     ax.set_title('Distribution of B by Chromosome', fontsize=15, fontweight='bold')
     ax.tick_params(axis='x', rotation=45, labelsize=10)
     ax.tick_params(axis='y', labelsize=10)

     plt.tight_layout()
     plt.savefig(output_path, dpi=300)
     print(f"Plot saved to {output_path}")
