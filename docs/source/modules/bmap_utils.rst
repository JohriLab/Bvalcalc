B-map Utilities for VCF
=========================

.. code-block:: bash

    Bvalcalc --Bmap Bmap.csv --positions variants.vcf

**-\-Bmap [path/to/Bmap.csv]**
  B-map utilities for getting B statistics for specific sites in a VCF/txt file. Provide path to B map as following argument.

Core Arguments
---------------

**-\-positions**
  VCF or two column CSV file [CHR,POS] with list of positions (typically variants) to get B for from the B-map.

Optional Arguments
-------------------

**-/-plot_distribution [path]**
  Output path for a plot of the distribution of B across each chromosome (default: `B_distribution.png`).
  
**-/-out**
**-/-out_minimum**
**-/-out_maximum**
**-/-quiet**


def parseBmapArgs(argv=None):
    parser = argparse.ArgumentParser(description="B-map utilities for getting B statistics for sites in a VCF/txt file with specific positions.")
    parser.add_argument('--positions', type=str, required=True, help="VCF or  input as following argument")
    parser.add_argument('--plot_distribution', nargs='?', const='B_distribution.png', default=None, help="")   
    parser.add_argument('--out', type=str, default=None,
                        help="Path to save per-site B for variant sites in the VCF/txt file. Results are not saved if --out is not specified.")
    parser.add_argument('--out_minimum', action='store_true', help="If set, only the sites ABOVE the given threshold of B will be returned, i.e. B > [threshold].")   
    parser.add_argument('--out_maximum', action='store_true', help="If set, only the sites BELOW the given threshold of B will be returned, i.e. B < [threshold].")   
    parser.add_argument('--quiet', action='store_true', help="If set, silence print statements.")
    return parser.parse_args(argv)