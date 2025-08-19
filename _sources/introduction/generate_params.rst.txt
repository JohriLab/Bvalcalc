Generate Parameters
=================================

Population parameter files can be retrieved from one of the pre-built templates with the following command. The templates are example parameters for your convenience and should be tailored to your population of interest where appropriate. 

To learn more about what the parameters represent and how to optimize them for your analysis see the guide :doc:`Tailoring Parameters <../guides/params>`.

Arguments
----------

**-\-generate_params [SPECIES]**
    Save a local copy of a population parameter file from a chosen species template.
    
    Available species: ``human``, ``drosophila``, ``arabidopsis``, ``mouse``, ``pfalciparum``, ``celegans``, ``selfing``, ``dromel_cds``, ``dromel_utr``, ``dromel_phastcons``

**-\-dir [path/to/NewParams.py]**
    Optional path to save parameters (defaults to `./SpeciesParams.py`).

Example
--------

.. code-block:: bash

    Bvalcalc --generate_params human

This will save a copy of `HumanParams.py` with population parameters specific to humans in the current folder, ready for calculating B
