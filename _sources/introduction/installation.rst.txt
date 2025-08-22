Installation
============


Install **Bvalcalc** in a new conda environment via conda-forge (recommended):

.. code-block:: bash

   conda create -n bvalcalc -c conda-forge bvalcalc
   conda activate bvalcalc

Verify the install worked by running the CLI "help":

.. code-block:: bash

   Bvalcalc --help

You should see usage modes listed with e.g. ``Welcome to Bvalcalc v0.6.5!``

.. note::
   Bvalcalc supports both ``bvalcalc`` (lowercase) and ``Bvalcalc`` (capitalized) for the CLI.
   For the Python API, use uppercase: ``import Bvalcalc``.

PyPI
------------

Install **Bvalcalc** via pip from PyPI:

.. code-block:: bash

   pip install bvalcalc



If you encounter dependency conflicts or prefer to install from PyPI with specific dependency versions, use this manual approach:

.. code-block:: bash

   conda create -n bvalcalc-env python=3.10 -y
   conda activate bvalcalc-env
   conda install -c conda-forge "numpy>=1.26" "scipy>=1.13" matplotlib -y
   pip install bvalcalc

Dependencies
------------

- Python ≥ 3.10
- NumPy  
- SciPy  

These should be pulled in automatically, but you can upgrade them manually if needed:

.. code-block:: bash

   pip install --upgrade numpy scipy
