Calculate Gene B Recovery
=========================

.. code-block:: bash

    Bvalcalc --gene --params Params.py

**-\-gene**
  Calculate B values for a range of neutral sites flanking a single selected element.

Core Arguments
--------------

**-\-params [path/to/YourParams.py]** 
  Path to a Python file defining population genetic parameters, see here for accessing pre-built templates, :doc:`Generate Parameters <../introduction/generate_params>`, and here for adjusting parameters to new species/populations, :doc:`Tailoring Parameters <../guides/params>`.

**-\-element_size [int]**
  Total length of the selected region, e.g. gene or exon, (default: `10000`).

**-\-flank_len [int]**
  Length of the flanking neutral region to calculate B (default: `40000`)

Optional Arguments
------------------

**-\-plot [path]**  
  Generate a B recovery slope output with a specified path (default: `./Bplot.png` if no path is given)

**-\-out [path]**  
  Write B-values to a CSV file with the specified path (must also provide ``--out_binsize``)

**-\-out_binsize [int]**  
  Bin size to average B-values in the CSV output, required if ``--out`` is used.

**-\-pop_change**
  If included, compute current B (``Bcur``) under a step population size change, as described in `Johri et al. (2021) <https://doi.org/10.1093/molbev/msab050>`_. 
  Note that ``Bcur`` and ``time_of_change`` should be set in the parameters file when active.

**-\-gamma_dfe**
  If included, use a gamma distribution to define the DFE (instead of fixed ``f0``, ``f1``, ``f2``, ``f3``). 
  Note that ``mean``, ``shape`` and ``proportion_synonymous`` should be set in the parameters file when active.

**-\-constant_dfe**
  If included, use a constant fixed ``s`` value as the DFE of selected sites (instead of fixed ``f0``, ``f1``, ``f2``, ``f3``). 
  Note that ``s`` and ``proportion_synonymous`` should be set in the parameters file when active.
  
**-\-quiet**
  Suppress console output

Example
-------

.. code-block:: bash

    Bvalcalc --gene \
        --params HomSap_Cds_Params.py \
        --element_size 10000 \
        --flank_len 40000 \
        --plot ./Bplot.png

    # B for adjacent site: 0.9617168921951618
    # Mean B for flanking region: 0.979953881508711

Calculates B for a 40 kb neutral region flanking a single gene under selection of length 10 kb using template human CDS parameters, and plot the output.

.. image:: /_static/images/gene_Bplot.png
   :alt: B region example
   :class: with-shadow
   :align: center