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

**-\-plot_distribution [path]**
  Output path for a plot of the distribution of B across each chromosome (default: `B_distribution.png`).
  
**-\-out [path]**
  Path to save per-site B for variant sites in the VCF/txt file. Results are not saved if ``--out`` is not specified.

**-\-out_minimum [float]**
  Filters positions by B so that only the sites ABOVE the given threshold of B will be returned, i.e. when B >= [out_minimum] is True.
  
**-\-out_maximum [float]**
  Filters positions by B so that only the sites BELOW the given threshold of B will be returned, i.e. when B <= [out_maximum] is True.

**-\-quiet**
  Silence print statements.

Example
-----------

.. code-block:: bash

    Bvalcalc --Bmap Bmap.csv \
      --positions variants.vcf \
      --out variants_B_above_0.9.csv \
      --out_minimum 0.9

This will save a csv file listing the position and B of each site in the VCF with B above or equal 0.9. These variants represent the most neutrally-evolving alleles in the VCF and may be useful for demographic inference, see :doc:`Demographic Inference with B-map <../guides/demography>`.