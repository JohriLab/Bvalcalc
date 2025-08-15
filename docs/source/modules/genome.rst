Calculate Genome B-map
===============================

.. code-block:: bash

    Bvalcalc --genome --pop_params Params.py --bedgff_path CDS.bed

**-\-genome**
    Calculate a B-map for all neutral and conserved sites across the genome, considering linked and unlinked effects of selection from all conserved elements.

Core Arguments
--------------

**-\-pop_params [path/to/YourParams.py]** 
  Path to a Python file defining population genetic parameters, see here for accessing pre-built templates, :doc:`Generate Parameters <../introduction/generate_params>`, and here for adjusting parameters to new species/populations, :doc:`Tailoring Parameters <../guides/params>`.

**-\-bedgff_path [path/to/example.bed]**  
    Path to an annotation file of selected elements, in BED, GFF3 or CSV format (CHR,START,END)

Recommended Arguments
---------------------

**-\-chr_sizes [path/to/chr_sizes.csv]**  
    Path to a file specifying chromosome sizes in CSV format (CHR,END). If not provided, chromosome boundaries will default to the end of the last annotated gene

**-\-rec_map [path/to/rec_map.csv]**  
    Optional recombination (crossover) map in CSV format (CHR,START,RATE), where "RATE" is a multiplication factor for ``r`` in the parameters file. Note that recombination rates are averaged over chunks.

**-\-gc_map [path/to/gc_map.csv]**  
    Optional gene conversion initiation rate map in CSV format (CHR,START,RATE), where "RATE" is a multiplication factor for ``g`` in the parameters file. Note that the same map can be used for crossover and gene conversion rates, rates are averaged over chunks.

**-\-pop_change**
    If included, compute current B under a step population size change, as described in `Johri et al. (2021) <https://doi.org/10.1093/molbev/msab050>`_. Note that ``Ncur`` and ``time_of_change`` should be set in the parameters file when active.

Optional Arguments
------------------

**-\-out [path]**  
    Write B-values to a CSV file with the specified path (must also provide ``--out_binsize``)

**-\-out_binsize [int]**  
    Bin size to average B-values in the CSV output, required if ``--out`` is used.

**-\-gamma_dfe**
  If included, use a gamma distribution to define the DFE (instead of fixed ``f0``, ``f1``, ``f2``, ``f3``). 
  Note that ``mean``, ``shape`` and ``proportion_synonymous`` should be set in the parameters file when active.

**-\-constant_dfe**
  If included, use a constant fixed ``s`` value as the DFE of selected sites (instead of fixed ``f0``, ``f1``, ``f2``, ``f3``). 
  Note that ``s`` and ``proportion_synonymous`` should be set in the parameters file when active.
  
**-\-hri**
    If included, will enable post-hoc calculation of B under HRI (B'; Becher and Charlesworth 2025), for low recombination regions. By default, classic B values are used in these regions.

**-\-prior_Bmap [path/to/prior_Bmap.csv]**  
    Optional prior B-value map (``.csv`` format). Used to multiply the newly calculated B-values by a per-site prior (e.g. for regions under different selection parameters). Format: ``Chromosome,Start,Conserved,B``. Note that ``Conserved`` is required for parsing but does not affect output

**-\-chunk_size [int]**  
    Size of chunks to process in each B calculation (default: ``20000``). It may be useful to increase this in large chromosomes with sparse selection for tractability, though consider how this may affect analysis in conjunction with the number of ``precise_chunks``.

**-\-precise_chunks [int]**  
    Number of chunks on either side of a focal chunk to calculate precisely (default: ``3``). Increasing this beyond the default will lead to more precise results though reduces tractability, consider how this may affect analysis in conjunction with the ``chunk_size``.

**-\-verbose**  
    Print per-chunk processing summaries (default: `False`)

**-\-quiet**  
    Suppress console output

Example
-------

.. code-block:: bash

    Bvalcalc --genome \
      --pop_params DrosophilaParams.py \
      --bedgff_path drosophila_CDS.bed \
      --chr_sizes chrom_sizes.txt \
      --out Bvalues_drosophila_CDS.csv \
      --out_binsize 1000

Calculates a B-map for across the genome considering all CDS regions. Output of B values in 1kb bins for the region will be saved.

Notes
------

A caveat to the ``--region`` and ``--genome`` modes is that by default they combine and simplify distant elements in discrete chunks which can slightly change the distance of distant conserved elements when
calculating B. The default chunk size is 20kb and the window within which calculations are perfectly precise is three chunks in each direction (140kb total). This allows for vastly improved performance
and typically will not result in directional biases of B estimates for most analyses. 

To achieve more exact results you can specify the size of the chunks with ``--chunk_size``, and the size of the window to
perform perfectly precise calculations with ``--precise_chunks``, though this will come at the cost of perfomance so consider using HPC resources or limiting to a specific region with ``--region``.
