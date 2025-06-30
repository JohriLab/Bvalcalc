Demographic Inference with B-map
========================================

Demographic inference is confounded by the linked effects of selection; see `Johri et al. (2021) <https://doi.org/10.1093/molbev/msab050>`_ for biases in MSMC-like and SFS-based approaches, and `Marsh and Johri. (2024) <https://doi.org/10.1093/molbev/msae118>`_ for biases in ARG-based approaches.

To best avoid biases, it is essential to only use the most neutrally-evolving sites for demographic inference , i.e. the sites least affected by sweeps and BGS (highest B). 

Using ``--Bmap`` (see :doc:`B-map Utilities for VCF <../modules/bmap_utils>`), we can use a B-map from ``--genome`` to filter a VCF or CSV to keep only the most neutrally evolving sites.

Plotting the B distribution
----------------------------

The following command will find the B of each position in the VCF using the B-map, print a brief summary to standard out and save a simple plot (`B_distribution.png`) of the results.

.. code-block:: bash

    Bvalcalc --Bmap your_B_map.csv \ 
        --positions your.vcf \
        --plot_distribution
        --out variants_B.csv

Now we can open up `B_distribution.png` which will help decide on a minimum B cut-off for our demographic inference analysis (e.g. B >= 0.9).

In addition we have `variants_B.csv` which lists the B for each position in case you'd like to run your own stats or do your own plotting of the results.

Saving positions with high B
-----------------------------

Next, we can pick a cut-off and filter our list of positions to keep only sites with B >= 0.9 and save the output to `filtered_positions.csv`.

.. code-block:: bash

    Bvalcalc --Bmap your_B_map.csv \ 
        --positions your.vcf \
        --out_minimum 0.9
        --out filtered_positions.txt
        --bcftools_format

Adding the ``--bcftools_format`` option removes the B value column and reformats the output so it's easier to filter with bcftools.

Filtering a VCF
-----------------

Using the newly saved `filtered_positions.csv` and our original VCF we can filter using bcftools:

.. code-block:: bash

    bcftools view \
        -R variants_B_above_0.9.txt \
        your.vcf -Ov -o filtered.vcf

Now the `filtered.vcf` only contains sites with B >= 0.9 and is ready for more accurate demographic inference with your tool of choice!