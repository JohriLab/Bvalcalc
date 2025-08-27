Generate Parameters
=================================

Population parameter files can be retrieved from one of the pre-built templates with the following command. The templates are example parameters for your convenience and should be tailored to your population of interest where appropriate. 

To learn more about what the parameters represent and how to optimize them for your analysis see the guide :doc:`Tailoring Parameters <../guides/params>`.

Arguments
----------

**-\-generate_params [SPECIES]**
    Save a local copy of a population parameter file from a chosen species template in format ``GenSpe_Elementtype_Params.py``.
    
    Available species: ``selfing``, ``aratha_cds``, ``aratha_phastcons``, ``dromel_cds``, ``dromel_utr``, ``dromel_phastcons``, ``homsap_cds``

**-\-dir [path/to/NewParams.py]**
    Optional path to directory to save parameters in (defaults to current directory).

Example
--------

.. code-block:: bash

    Bvalcalc --generate_params homsap_cds

This will save a copy of `HomSap_Cds_Params.py` with population parameters specific to human CDS (coding sequence) in the current folder, ready for calculating B
