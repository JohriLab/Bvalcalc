Quickstart Tutorial
====================

This tutorial (5-10 minutes) walks through the core functionality of calculating B with Bvalcalc!

Installation
------------

Bvalcalc is a python program called from the command line (terminal, shell, zsh, bash etc.)

You can install Bvalcalc on the command line via pip:

.. code-block:: bash

   pip install Bvalcalc

Specifying popgen parameters
------------------------------

To calculate B, we need evolutionary information about the species of the population in a parameters file.

Copy population genetic parameters from one of the templates, e.g. Drosophila.  
In your own analysis you'll need to check the literature and use informed values for your population.

.. code-block:: bash

   Bvalcalc --generate_params drosophila

Open ``./DrosophilaParams.py`` in your text editor of choice. This file contains example popgen parameters that allow us to accurately calculate B.

Calculating a B value
----------------------

Let's start by calculating B at a site 1kb from the edge of a single conserved element 5kb in length.

.. code-block:: bash

   Bvalcalc \
   --site \
   --pop_params ./DrosophilaParams.py \
   --gene_size 5000 \
   --distance 1000

The B-value for this site should be printed in the console as just above 0.99, which indicates diversity is expecteed to be reduced 1% due to BGS.

Calculating B recovery from a single element
---------------------------------------------

Now, let's calculate the recovery of B as a function of distance from the 5kb element so we can plot it (gene_B.png).

.. code-block:: bash

   Bvalcalc \
   --gene \
   --gene_size 5000 \
   --pop_params ./DrosophilaParams.py \
   --plot_output gene_B.png

Have a look at the plot and the results printed to the console, you'll notice B decays with distance from the selected element. It's still a modest reduction (<3.1%), but remember, across a genome ALL selected elements will contribute to B at any given site.

Calculating B for a region of the genome
-----------------------------------------

Alright, now let's look at part of a chromosome.  
We can use a BED file (or GFF/CSV) that specifies which genomic ranges are conserved to calculate B for a region in the genome.  
Let's calculate B for a 1Mb region in the middle of chromosome 2R [9500000-10500000] and save it to 1Mb_B.png.

.. code-block:: bash

   Bvalcalc \
   --region 2R:9500000-10500000 \
   --pop_params ./DrosophilaParams.py \
   --bedgff_path examples/dmel6_2R_genes.csv \
   --plot_output 1Mb_B.png

Have a look at the plots: the blue sections of the graph indicate neutral regions and black indicates conserved elements.  
That's all that's necessary for many analyses, especially if you're only interested in B values for a specific region of the genome, or are testing against simulated results.

Calculating a complete B-map
-----------------------------

If you wanted to generate a complete B-map for all sites across all chromosomes you would use the following command, though note it's a lot more data to crunch and maps are already available for Drosophila so no need to run it!

.. code-block:: bash

   Bvalcalc \
   --genome \
   --pop_params ./DrosophilaParams.py \
   --bedgff_path examples/dmel6_2R_genes.csv \
   --out Dmel_Bmap.csv \
   --out_binsize 1000

If you had run that command, you'd get a B-map! 

B-maps are useful to identify highly conserved regions of the genome, as a null-model for inference, e.g. :doc:`SweepFinder2 with B-map <./sweepfinder2>`, or to select the most neutrally-evolving sites for e.g. demographic inference, see :doc:`Demographic Inferece with B-map <./demography>`.
