Calculate B
===============

The Python API can be used to access core B calculation functions that take vectorized inputs.

Python functions
-----------------

.. autofunction:: Bvalcalc.get_params

.. autofunction:: Bvalcalc.calculateB_linear

.. autofunction:: Bvalcalc.calculateB_unlinked

Example
---------

.. code-block:: python

   from Bvalcalc import get_params, calculateB_unlinked, calculateB_linear

   # Path to your Params.py file
   params = get_params("./DrosophilaParams.py")

   calculateB_unlinked(unlinked_L = 200000, params = params)
   # 0.9998476570803477

   calculateB_linear(distance_to_element = 500, length_of_element = 10000, params = params)
   # array(0.98721754)

This will import the relevant functions, get the popgen parameters from a relevant params file, see :doc:`Generate Parameters <../introduction/generate_params>`. 
Then B will be calculated from 200kb of unlinked sites using ``calculateB_unlinked``, and B from a linked selected element of length 10kb, 500bp away is calculated using ``calculateB_linear``.

Notes
---------

Note that ``calculateB_linear`` assumes a consistent crossover and gene conversion rate across both the length of and distance to the selected element, in the CLI, variable recombination rates are accounted for with 
the more complex function ``calculateB_recmap``; if accessing this function is important for your work feel free to raise it as an issue on GitHub as I could add it as a public API.