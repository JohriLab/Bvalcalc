Calculate Genome B-map
===============================

Calculate a B-map for all neutral and conserved sites across the genome, considering linked and unlinked effects of selection from all conserved elements.


.. code-block:: console

    bvalcalc --genome --pop_params YourParams.py --bedgff_path annotations.gff --chr_sizes chrom_sizes.txt

Core Arguments
--------------

**--pop_params**  
Path to a Python file defining population genetic parameters  
*(e.g. `Params.py`, see [Population Parameters](../your-reference-page))*

**--bedgff_path**  
Path to an annotation file of selected elements in BED or GFF3 format

**--chr_sizes**  
Path to a tab-delimited chromosome size file  
If not provided, chromosome boundaries will default to the end of the last annotated gene

Optional Arguments
------------------

**--chunk_size**  
Size of genomic chunks processed per thread (default: `20000`)

**--precise_chunks**  
Number of chunks on either side of a focal chunk to calculate precisely (default: `3`)

**--pop_change**  
Apply a step-change population size history (defined by `Ncur` and `time_of_change`) to output `Bcur` instead of `Banc`

**--gamma_dfe**  
Use a gamma-distributed DFE to replace fixed `f0–f3` proportions

**--prior_Bmap**  
Optional input CSV defining prior B-values per site (e.g. combining multiple annotation layers)  
Format: `Chromosome,Position,Conserved,B`  
→ `Conserved` is required for parsing but doesn't affect the result

**--rec_map**  
Optional recombination rate map in CSV format with header: `start,rate`  
Rates are averaged within each chunk

**--gc_map**  
Optional gene conversion rate map, same format as recombination map

**--neutral_only**  
If used, any plots generated will include only neutral (non-conserved) sites

**--out [path]**  
Write B-values to a `.csv` file (must also use `--out_binsize`)

**--out_binsize**  
Bin size for averaging B-values in CSV output. Required when `--out` is used

**--verbose**  
Print detailed progress and per-chunk summaries

**--quiet**  
Suppress console output

Example
-------

.. code-block:: console

    bvalcalc --genome \
      --pop_params HumanParams.py \
      --bedgff_path annotations.gff \
      --chr_sizes chrom_sizes.txt \
      --out Bvalues_chr2R.csv \
      --out_binsize 1000

Notes
-----

- If `--out` is used, you **must** also provide `--out_binsize`
- B-values are calculated per base by default unless binning is specified
- Recombination and GC rates are averaged within each chunk boundary
- Use `--neutral_only` to visualize only unconserved sites in plots
