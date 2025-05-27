Calculate B
===============

The Python API can be used to access core B calculation functions that take vectorized inputs.

Python functions
-----------------

.. autofunction:: bvalcalc.get_params

.. autofunction:: bvalcalc.calculateB_linear

.. autofunction:: bvalcalc.calculateB_unlinked

Example
---------

.. code-block:: python

   from bvalcalc import get_params, calculateB_unlinked, calculateB_linear

   params = get_params("bvalcalc/templates/DrosophilaParams.py")

   calculateB_unlinked(unlinked_L = 200000, params = params)

   calculateB_linear(distance_to_element = 500, length_of_element = 10000, params = params)

This will import the relevant functions, get the popgen parameters from a relevant params file, see :doc:`Generate Parameters <../introduction/generate_params>`. 
Then B will be calculated from 200kb of unlinked sites using ``calculateB_unlinked``, and B from a linked selected element of length 10kb, 500bp away is calculated using ``calculateB_linear``.

Notes
---------

Note that ``calculateB_linear`` assumes a consistent crossover and gene conversion rate across both the length of and distance to the selected element, in CLI, variable recombination rates are accounted for with 
the more complex function ``calculateB_recmap``; if accessing this function is important for your work feel free to raise it as an issue on GitHub as I could add it as a public API.