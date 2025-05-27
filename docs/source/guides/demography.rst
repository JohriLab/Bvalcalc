Demographic Inference with B-map
========================================

Demographic inference tools are confounded by the linked effects of selection; see Johri et al. (2021) for biases in SFS-based approaches, and Marsh et al. (2024) for biases in ARG-based approaches.

To best avoid biases, it is essential to only use the most neutrally-evolving sites for demographic inference , i.e. the sites least affected by sweeps and BGS (highest B). 

Using the ``--Bmap`` utilities mode, we can get the B of a list of positions from a VCF or CSV from a B-map (e.g. from ``--genome``) to help select the most neutrally evolving sites.

Plotting the B Distribition
----------------------------

The following command will find the B of each position in the VCF using the B-map, print a brief summary to standard out and save a simple plot (`B_distribution.png`) of the results.

.. code-block:: bash

    Bvalcalc --Bmap your_B_map.csv \ 
        --positions your.vcf \
        --plot_distribution
        --out variants_B.csv

Now we can open up `B_distribution.png` which will help decide on a minimum B cut-off for our demographic inference analysis (e.g. B >= 0.9).

In addition we have `variants_B.csv` which lists the B for each position in case you'd like to run your own stats or do your own plotting of the results.

Saving Positions with high B
-----------------------------

Next, we can pick a cut-off and filter our list of positions to keep only sites with B >= 0.9 and save the output, we add the ``bcftools_format`` so it's easier to filter in the next step.

.. code-block:: bash

    Bvalcalc --Bmap your_B_map.csv \ 
        --positions your.vcf \
        --out_minimum 0.9
        --out filtered_positions.txt
        --bcftools_format

This will save `filtered_positions.csv` which lists positions with B >= 0.9 in the right format to filter our VCF with a tool like ``bcftools view``

Filtering a VCF
-----------------

Using the newly saved `filtered_positions.csv` and our original VCF:

.. code-block:: bash

    bcftools view \
        -R variants_B_above_0.9.txt \
        your.vcf -Ov -o filtered.vcf

Now your `filtered.vcf` only contains sites with B >= 0.9 and is ready for more accurate demographic inference with your tool of choice!