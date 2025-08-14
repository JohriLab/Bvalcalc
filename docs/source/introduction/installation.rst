Installation
============

PyPI
------------

Install **bvalcalc** via pip from PyPI:

.. code-block:: bash

   pip install bvalcalc

... or for the dev version!

.. code-block:: bash

   pip install --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple --upgrade bvalcalc

If it's not working try this:

.. code-block:: bash

   conda create -n bvalcalc-env python=3.10 -y
   conda activate bvalcalc-env
   conda install -c conda-forge "numpy>=1.26" "scipy>=1.13" matplotlib -y
   pip install --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple --upgrade bvalcalc

Verify the install worked by running the CLI "help":

.. code-block:: bash

   bvalcalc --help

You should see usage modes listed (``--site``, ``--gene``, etc.).

Dependencies
------------

- Python ≥ 3.10
- NumPy  
- SciPy  

These should be pulled in automatically, but you can upgrade them manually if needed:

.. code-block:: bash

   pip install --upgrade numpy scipy
