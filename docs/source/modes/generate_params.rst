Generating population parameters
=========================

A population parameter file can be generated from one of the pre-built templates in `Bvalcalc` using `--generate_params`.

This file provides population genetic parameters such as the effective population size and per-base mutation/recombination rates which are all variables that impact the linked effects of selection (i.e., B).

.. code-block:: console

    bvalcalc --generate_params human

This will save a copy of `SpeciesParams.py` with population parameters specific to humans in current folder (can specify with `--dir`).

Available templates: `human`, `drosophila`, `arabidopsis`, `mouse`, `pfalciparum`, `celegans`, `selfing`

These are example parameters for your convenience and should be tailored to your population of interest where appropriate.

Note that the `selfing`, `pfalciparum` and `arabidopsis` templates have an additional parameter for self-fertilization (F; Wright's inbreeding coefficient), which adjusts some values.
