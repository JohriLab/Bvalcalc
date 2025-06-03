Calculate Gene B Recovery
=========================

.. code-block:: bash

    Bvalcalc --gene --pop_params Params.py

**-\-gene**
  Calculate B values for a range of neutral sites flanking a single selected element.

Core Arguments
--------------

**-\-pop_params [path/to/YourParams.py]** 
  Path to a Python file defining population genetic parameters, see [here for generating a pre-built template] and [here for tailoring your own parameters]

**-\-gene_size [int]**
  Total length of the selected region, e.g. gene or exon, (default: `10000`)

**-\-flank_len [int]**
  Length of the flanking neutral region to calculate B (default: `40000`)

Optional Arguments
------------------

**-\-plot_output [path]**  
  Generate a B recovery slope output with a specified path (default: `./Bplot.png` if no path is given)

**-\-out [path]**  
  Write B-values to a CSV file with the specified path (must also provide ``--out_binsize``)

**-\-out_binsize [int]**  
  Bin size to average B-values in the CSV output, required if ``--out`` is used.

**-\-pop_change**
  If included, compute current B (``Bcur``) under a step population size change, see the associated guide on [demography]. 
  Note that ``Bcur`` and ``time_of_change`` should be set in the parameters file when active.

**-\-gamma_dfe**
  If included, use a gamma distribution to define the DFE (instead of fixed ``f0``, ``f1``, ``f2``, ``f3``). 
  Note that ``mean``, ``shape`` and ``proportion_synonymous`` should be set in the parameters file when active.

**-\-quiet**
  Suppress console output

Example
-------

.. code-block:: bash

    Bvalcalc --gene \
        --pop_params HumanParams.py \
        --gene_size 10000 \
        --flank_len 40000 \
        --plot_output ./Bplot.png

Calculates B for a 40 kb neutral region flanking a single gene under selection of length 10 kb using example human parameters, and plot the output.