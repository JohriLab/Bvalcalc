Calculate B
===============

The Python API can be used to access core B calculation functions that take vectorized inputs (numpy arrays). Here we show how to calculate B from linked sites (``calculateB_linear``) and unlinked sites (``calculateB_unlinked``).

Python functions
-----------------

.. autofunction:: Bvalcalc.get_params

.. autofunction:: Bvalcalc.calculateB_linear

.. autofunction:: Bvalcalc.calculateB_unlinked

.. autofunction:: Bvalcalc.calculateB_hri

Example
---------

.. code-block:: python

   from Bvalcalc import get_params, calculateB_unlinked, calculateB_linear, calculateB_hri

   # Path to your Params.py file
   params = get_params("./DroMel_Cds_Params.py")

   calculateB_unlinked(unlinked_L = 200000, params = params)
   # np.float64(0.9996953434292404)

   calculateB_linear(distance_to_element = 500, length_of_element = 10000, params = params)
   # array(0.98721754)

   calculateB_hri(distant_B = 0.7, interfering_L = 10000, params = params)
   # np.float64(0.22627504209854507)

This will import the relevant functions, get the popgen parameters from a relevant params file, see :doc:`Generate Parameters <../introduction/generate_params>`. 
Then B will be calculated from 200kb of unlinked sites using ``calculateB_unlinked``, and B from a linked selected element of length 10kb, 500bp away is calculated using ``calculateB_linear``.
Finally, B' is calculated for a non-recombining region containing 10000bp of selected sites, with a prior B of 0.7 from distant selected sites outside the region.

Notes
---------

Note that ``calculateB_linear`` assumes a consistent crossover and gene conversion rate across both the length of and distance to the selected element, in the CLI, variable recombination rates are accounted for with 
the more complex function ``calculateB_recmap``; if accessing this function is important for your work feel free to raise it as an issue on GitHub as I could describe its API usage.