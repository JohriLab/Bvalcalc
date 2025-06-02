SweepFinder2 with B-map
============================

**SweepFinder2** (**SF2**; DeGiorgio et al. 2016) is a population genetic program for performing genome-wide scans to identify selective sweeps, which uses a B-map to avoid confounding effects of BGS.

A B-map calculated by **Bvalcalc**'s ``--genome`` or ``--region`` module can be used as the B-value file input to **SF2** (see 3.3. B-value file in the **SF2** manual).

SF2 B-map format
---------------------

The required format for **SF2** is a tab-delimited file for a single chromosome with a header (position   bvalue), e.g.:

.. code-block:: console

    # SF2 B-value file format
    position    bvalue
    340001  0.2
    560010  0.5
    780210  0.95

    # Bvalcalc B-map format
    Chromosome,Position,B
    chr1,1,0.7
    chr1,1001,0.6
    chr1,2001,0.5

Note that this is different from the B-map format output of **Bvalcalc** which is a CSV and has an additional column with "chromosome". Also, it requires pulling out the B-values of specific positions, rather than a genome-wide map.

Bvalcalc to SF2 format
-------------------------

First, assuming you have a B-map from **Bvalcalc** you'll need to get the B-values for your variant positions (the SF2 grid file).

.. code-block:: bash

    # Add a chromosome column to your grid file and save to a new file
    sed 's/^/chr1,/' chr1_grid.txt > chr1_positions.csv

    # Then get their B-values from the B-map
    Bvalcalc --Bmap B_map.csv \
        --positions chr1_positions.csv \
        --out chr1_bvalues.csv

This will save the B-values for your list of chromosome 1 positions (chr1_positions.csv) to chr1_bvalues.csv. 

Now, using awk we can pull out positions for a single chromosome and reformat it to the SF2 B-value format

.. code-block:: bash

    # Add a header to a new file
    echo "position\tbvalue" > chr1_bvalues.tsv
    
    # Reformat for SF2-format and append
    awk -F, 'BEGIN { OFS="\t" } \
     $1=="chr1" { print $2,$3 }' chr1_bvalues.csv >> chr1_bvalues.tsv

Now you have chr1_bvalues.tsv ready as input to SF2!

Notes
------

**SF2** analyses typically require a recombination map, so you can use the same recombination map when calculating the B-map by adding ``--rec_map your.map``, just note you'll need to convert between centiMorgan and rate format.


