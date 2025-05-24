Generate Parameters
=================================

Population parameter files can be retrieved from one of the pre-built templates with the following commands.

**-\-generate_params [SPECIES]**
    Save a local copy of a population parameter file from a chosen species template.
    
    Available species: `human`, `drosophila`, `arabidopsis`, `mouse`, `pfalciparum`, `celegans`, `selfing`

**-\-dir [path/to/NewParams.py]**
    Optional path to save parameters (defaults to `./SpeciesParams.py`).

The python file provides parameters that are used to calculate B such as the effective population size and per-base mutation/recombination rates.

These are example parameters for your convenience and should be tailored to your population of interest where appropriate.

To learn more about what the parameters represent and how to optimize them for your analysis see the guide see :doc:`Tailoring Parameters <../guides/params>`.

Example
--------

.. code-block:: bash

    Bvalcalc --generate_params human

This will save a copy of `HumanParams.py` with population parameters specific to humans in the current folder
