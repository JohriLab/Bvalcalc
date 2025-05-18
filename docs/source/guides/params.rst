Tailoring Parameters
=========================


Selfing species
---------------

Self-fertilizing (selfing) species have different evolutionary dynamics from obligate outcrosses which impact BGS (i.e., B).

The parameters file can be modified for selfing populations by altering population size (Nanc/Ncur), crossover rate (r), gene conversion rate (g) and dominance coefficient by an additional parameter: f, which is Wright's inbreeding coefficient (F). F can be calculated from the selfing rate (F = S/(2-S)). Note that the `arabidopsis` and `pfalciparum` default templates have the inbreeding parameter included.

For analysis of selfing populations, we recommend tailoring parameters from the `selfing` template.

.. code-block:: bash

    Bvalcalc --generate_params selfing