BED/GFF Input
=============

**Bvalcalc** estimates the effects of BGS on linked (and unlinked) sites due to direct purifying selection in conserved regions of the genome.

When using the ``--region`` and ``genome`` modules, these conserved regions are indicated by a BED/GFF/CSV annotation file:

**-\-bedgff [path/to/example.bed]**  
    Path to an annotation file of selected elements, in BED, GFF3 or CSV format (CHR,START,END). Header information, indicated by ``#``, is ignored.

Annotations
------------

The annotations indicated in the BED/GFF/CSV input will all be assumed to have the same distribution of fitness effects (DFE), though this may not always be the case.

Typical annotations that are expected to experience purifying selection are exons and regulatory regions. Inferred conserved elements such as phastCons or GERP elements can make for excellent inputs for **Bvalcalc**.

If you have DFE information for the different annotation types (e.g., CDS and UTRs), consider running **Bvalcalc** separately for each and combining the results, see :doc:`Multiple DFEs <../guides/multiple_dfes>`.

Example
------------

.. code-block:: console

    #GFF3 format
    2R	.	.	326768	705848	.	.	.	.
    2R	.	.	854814	855571	.	.	.	.
    2R	.	.	924636	951252	.	.	.	.

    #BED format
    2R	326768	705848
    2R	854814	855571
    2R	924636	951252

    #CSV format
    2R,326768,705848
    2R,854814,855571
    2R,924636,951252

Three headerless example inputs to ``--bedgff`` that will be processed identically: