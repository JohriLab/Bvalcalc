Calculate Single Site B
=========================

.. code-block:: bash

    Bvalcalc --site --params Params.py

**-\-site**
  Calculate B for a single site at a specified distance from a single selected element

Core Arguments
------------------

**-\-params [path/to/YourParams.py]** 
  Path to a Python file defining population genetic parameters, see here for accessing pre-built templates, :doc:`Generate Parameters <../introduction/generate_params>`, and here for adjusting parameters to new species/populations, :doc:`Tailoring Parameters <../guides/params>`.

**-\-element_size [int]**
  Total length of the selected region, e.g. gene or CDS, (default: `10000`)

**-\-distance [int]**
  Distance from the focal neutral site to the edge of the selected element in bp (default: `1`). Note that this assumes recombinant distance increases linearly with physical distance.

Optional Arguments
------------------

**-\-pop_change**
  If included, compute current B (``Bcur``) under a step population size change, as described in `Johri et al. (2021) <https://doi.org/10.1093/molbev/msab050>`_. 
  Note that ``Bcur`` and ``time_of_change`` should be set in the parameters file when active.

**-\-gamma_dfe**
  If included, use a gamma distribution to define the DFE (instead of fixed ``f0``, ``f1``, ``f2``, ``f3``). 
  Note that ``mean``, ``shape`` and ``proportion_synonymous`` should be set in the parameters file when active.

**-\-constant_dfe**
  If included, use a constant fixed ``s`` value as the DFE of selected sites (instead of fixed ``f0``, ``f1``, ``f2``, ``f3``). 
  Note that ``s`` and ``proportion_synonymous`` should be set in the parameters file when active.

**-\-quiet**
  Suppress console output

Example
-------
.. code-block:: bash

    Bvalcalc --site --params HomSap_Cds_Params.py --distance 1500 --element_size 10000

    # B for site 1500bp away from 10000bp region: 0.964573141826118

Calculates B for a single site 1500bp away from a gene under selection of length 10kb using the template human CDS parameters.