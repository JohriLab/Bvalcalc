Generate Parameters
=================================

Population parameter files can be retrieved from one of the pre-built templates using `Bvalcalc --generate_params`.

The python file provides parameters such as the effective population size and per-base mutation/recombination rates which are all variables that impact the linked effects of selection (i.e., B).

--generate_params [SPECIES]
    Save a local copy of a population parameter file from a chosen species template.
    Available species: `human`, `drosophila`, `arabidopsis`, `mouse`, `pfalciparum`, `celegans`, `selfing`

--dir [path/to/NewParams.py]
    Optional path to save parameters (defaults to `./SpeciesParams.py`).

These are example parameters for your convenience and should be tailored to your population of interest where appropriate.

To learn more about what the parameters represent and how to optimize them for your analysis see the guide [Tailoring Parameters].

Example
--------

.. code-block:: bash

    Bvalcalc --generate_params human

This will save a copy of `HumanParams.py` with population parameters specific to humans in the current folder
