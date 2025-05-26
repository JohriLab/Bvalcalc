Multiple DFEs
===============

Different annotated genomic elements that are expected to be under evolutionary constraint (e.g. CDS, UTRs, promoters, enhancers etc.) will have different distributions of fitness of effects (DFEs).

To account for different DFEs that may affect different annotated regions, you can run **Bvalcalc** for each DFE iteratively with the help of the ``--prior_Bmap`` flag in the ``--genome`` mode.

**-\-prior_Bmap [path/to/prior_Bmap.csv]**  
    Optional prior B-value map (`.csv` format). Used to multiply the newly calculated B-values by a per-site prior (e.g. for regions under different selection parameters). 
    
    Format: `Chromosome,Position,Conserved,B`. Note that `Conserved` is required for parsing but does not affect output

Walkthrough
------------

Let's briefly walk through how you would generate a B-map that considers the combined effects of CDS (strongly deleterious DFE) and UTRs (moderately deleterious DFE).
