Installation
============

PyPI
------------

Install **Bvalcalc** via pip from PyPI:

.. code-block:: bash

   pip install bvalcalc

(or: ``python -m pip install bvalcalc``)

Verify the install worked by running the CLI “help”:

.. code-block:: bash

   bvalcalc --help

You should see usage modes listed (“--site”, “--gene”, etc.).

Dependencies
------------

- Python ≥ 3.10
- NumPy  
- SciPy  

These should be pulled in automatically, but you can upgrade them manually if needed:

.. code-block:: bash

   pip install --upgrade numpy scipy
