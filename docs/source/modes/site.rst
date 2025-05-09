Single Site Calculation
================

Calculate B for a single site at a specified distance from a single selected element:

.. code-block:: console

    bvalcalc --site --pop_params YourParams.py --distance 100 --gene_size 5000

Core Arguments
------------------

**--pop_params**: Path to a Python file defining population genetic parameters  
  *(e.g. `Params.py`, see [Population Parameters](../your-reference-page))*

**--distance**: Distance from the focal neutral site to the edge of the selected element (default: `1`)

**--gene_size**: Total length of the selected region (default: `10000`)

Optional Arguments
------------------

**--pop_change**: If included, compute current B (`Bcur`) under a step population size change

**--gamma_dfe**: If included, use a gamma distribution to define the DFE (instead of fixed `f0,f1,f2,f3`)

**--quiet**: Suppress console output

Example
-------

.. code-block:: console

    bvalcalc --site --pop_params HumanParams.py --distance 1500 --gene_size 10000
