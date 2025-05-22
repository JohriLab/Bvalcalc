Calculate Single Site B
=========================

.. code-block:: bash

    Bvalcalc --site --pop_params Params.py

**-\-site**
  Calculate B for a single site at a specified distance from a single selected element

Core Arguments
------------------

**-\-pop_params [path/to/YourParams.py]** 
  Path to a Python file defining population genetic parameters, see [here for generating a pre-built template] and [here for tailoring your own parameters]

**-\-gene_size [int]**
  Total length of the selected region, e.g. gene or CDS, (default: `10000`)

**-\-distance [int]**
  Distance from the focal neutral site to the edge of the selected element in bp (default: `1`). Note that this assumes recombinant distance increases linearly with physical distance.

Optional Arguments
------------------

**-\-pop_change**
  If included, compute current B (`Bcur`) under a step population size change, see the associated guide on [demography]. 
  Note that `Bcur` and `time_of_change` should be set in the parameters file when active.

**-\-gamma_dfe**
  If included, use a gamma distribution to define the DFE (instead of fixed `f0,f1,f2,f3`). 
  Note that `mean`, `shape` and `proportion_synonymous` should be set in the parameters file when active.

**-\-quiet**
  Suppress console output

Example
-------
.. code-block:: bash

    Bvalcalc --site --pop_params HumanParams.py --distance 1500 --gene_size 10000

Calculates B for a single site 1500bp away from a gene under selection of length 10kb using example human parameters.