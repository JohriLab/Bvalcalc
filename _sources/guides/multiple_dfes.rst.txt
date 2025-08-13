Multiple DFEs
===============

Different annotated genomic elements that are expected to be under evolutionary constraint (e.g. CDS, UTRs, promoters, enhancers etc.) will have different distributions of fitness of effects (DFEs).

To account for different DFEs that may affect different annotated regions, you can run **Bvalcalc** for each DFE iteratively with the help of the ``--prior_Bmap`` flag in the ``--genome`` mode.

Arguments
----------

**-\-prior_Bmap [path/to/prior_Bmap.csv]**  
    Optional prior B-value map (`.csv` format). Used to multiply the newly calculated B-values by a per-site prior (e.g. for regions under different selection parameters). 
    
    Format: `Chromosome,Start,B`.

Example Walkthrough
--------------------

Let's briefly walk through how you would generate a B-map that considers the combined effects of CDS (strongly deleterious DFE) and UTRs (moderately deleterious DFE).

1. Save a local a template params file using ``--generate_params`` for each of your annotations (CDS and UTR). See :doc:`Generate Parameters <../introduction/generate_params>`

.. code-block:: bash

    Bvalcalc --generate_params drosophila
    scp DrosophilaParams.py CdsParams.py
    scp CdsParams.py UtrParams.py

2. Tailor the CdsParams.py/UtrParams.py files with the respective DFEs. See :doc:`Tailoring Parameters <./params>`

\

3. Run ``--genome`` to create a B-map considering only the first annotation type (e.g. CDS). See :doc:`Calculate Genome B-map <../modules/genome>`

.. code-block:: bash

    Bvalcalc --genome \
        --pop_params CdsParams.py \
        --bedgff_path drosophila_CDS.bed \
        --out CDS_B.csv \
        --out_binsize 1000

4. Run ``--genome`` and pass in the first annotations (CDS) Bmap ``--prior_Bmap`` as an additional argument.

.. code-block:: bash

    Bvalcalc --genome \
        --pop_params UtrParams.py \
        --bedgff_path drosophila_UTR.bed \
        --prior_Bmap CDS_B.csv \
        --out UTR_CDS_B.csv \
        --out_binsize 1000

We now have a B-map (UTR_CDS_B.csv) that considers the effects of both CDS and UTRs.