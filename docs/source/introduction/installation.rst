Installation
============

PyPI
------------

Install **Bvalcalc** via pip from PyPI:

.. code-block:: bash

   pip install Bvalcalc

... or for the dev version!

.. code-block:: bash

   pip install --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple --upgrade Bvalcalc

Verify the install worked by running the CLI “help”:

.. code-block:: bash

   Bvalcalc --help

You should see usage modes listed (``--site``, ``--gene``, etc.).

Dependencies
------------

- Python ≥ 3.10
- NumPy  
- SciPy  

These should be pulled in automatically, but you can upgrade them manually if needed:

.. code-block:: bash

   pip install --upgrade numpy scipy
