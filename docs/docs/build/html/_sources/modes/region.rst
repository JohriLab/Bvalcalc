Calculate Region B-map
============================

Calculate a B-map for a specified chromosomal region, considering linked and unlinked effects of selection from all conserved elements across the genome.


.. code-block:: console

    bvalcalc --region 2R:9260000-11700000 --pop_params YourParams.py --bedgff_path annotations.gff 

Core Arguments
--------------

**--region**  
Chromosome and coordinate range to calculate B over  
Format: `CHR:START-END` (e.g. `2R:9260000-11700000`)  
If not set, defaults to the entire chromosome in the annotation file

**--pop_params**  
Path to a Python file defining population genetic parameters  
*(e.g. `Params.py`, see [Population Parameters](../your-reference-page))*

**--bedgff_path**  
Path to an annotation file of selected elements, in BED or GFF3 format

Optional Arguments
------------------

**--chunk_size**  
Size of chunks to process in each B calculation (default: `20000`)

**--precise_chunks**  
Number of chunks on either side of a focal chunk to calculate precisely (default: `3`)

**--pop_change**  
Account for a step-change population size history (`Nanc → Ncur`), computing `Bcur` rather than `Banc`

**--gamma_dfe**  
Use gamma-distributed DFE parameters instead of fixed `f0–f3`

**--prior_Bmap**  
Optional prior B-value map (`.csv` format). Used to multiply the newly calculated B-values by a per-site prior (e.g. for combined annotations)  
Format: `Chromosome,Position,Conserved,B`  
→ `Conserved` is required for parsing but does not affect output

**--rec_map**  
Optional recombination map. CSV with columns: `start,rate`

**--gc_map**  
Optional gene conversion map. CSV with columns: `start,rate`

**--plot_output [path]**  
Output path for a genome-wide or region B-value plot (default: `genome_plot.png`)

**--neutral_only**  
If used with `--plot_output`, only neutral sites will be shown in the plot

**--out [path]**  
Path to output `.csv` file containing B-values. If omitted, no file will be saved  
→ If `--out` is provided, `--out_binsize` **must** also be specified

**--out_binsize**  
Bin size used to average B-values for CSV output. Required if using `--out`

**--verbose**  
Print per-chunk processing summaries (default: False)

**--quiet**  
Suppress console output

Example
-------

.. code-block:: console

    bvalcalc --region 2R:9260000-11700000 \
      --pop_params HumanParams.py \
      --bedgff_path annotations.gff \
      --chunk_size 20000 \
      --plot_output BregionPlot.png

Notes
-----

- You **must** provide a valid region or ensure your BED/GFF file covers only the desired region.
- If you use `--out`, you **must** also include `--out_binsize`.
- Recombination and GC maps will be averaged over chunks.
