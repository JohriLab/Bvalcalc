Gene Recovery Calculation
=====================

Compute B-values for a range of neutral sites flanking a single selected element (e.g., a gene under selection).

.. code-block:: console

    bvalcalc --gene --pop_params YourParams.py --gene_size 10000 --flank_len 40000

Core Arguments
--------------

**--pop_params**  
Path to a Python file defining population genetic parameters  
*(e.g. `Params.py`, see [Population Parameters](../your-reference-page))*

**--gene_size**  
Length of the region (e.g. gene) under selection (default: `10000`)

**--flank_len**  
Length of the flanking neutral region to calculate B (default: `40000`)

Optional Arguments
------------------

**--pop_change**  
Account for a step-change population size history, detailed by **Ncur** and **time_of_change** in the params file

**--gamma_dfe**
If included, use a gamma distribution to define the DFE (instead of fixed `f0,f1,f2,f3`)

**--plot_output [path]**  
Generate a B recovery slope output with a specified path (default: `./Bplot.png` if no path is given)

**--out [path]**  
Write B-values to a CSV file with the specified path (must also provide `--out_binsize`)

**--out_binsize**  
Bin size to average B-values in the CSV output, required if `--out` is used.

**--quiet**  
Suppress console output

Example
-------

.. code-block:: console

    bvalcalc --gene --pop_params HumanParams.py --gene_size 10000 --flank_len 40000 --plot_output ./Bplot.png
