Tailoring Parameters
=========================

When using **Bvalcalc** to study populations beyond the pre-built templates, it's important to tailor the popgen parameters to your population of interest.

To save a local copy of one of the templates, see :doc:`Generate Parameters <../introduction/generate_params>` and open up the Params.py file in your text editor of choice.

Core parameters
----------------
``x`` 
    Scaling factor that modifies ``N``, ``u``, ``r``, and ``g``. Keep as 1 for empirical analysis, only relevant for calculating B to compare against rescaled simulations. See `Marsh, Kaushik and Johri 2025 <https://doi.org/10.1101/2025.04.24.650500>`_.
``Nanc``
    Ancestral population size, often reported in literature. This is the population size value that scales the strength of selection. Can be roughly estimated from nucleotide diversity at neutrally evolving sites, given a mutation rate and no demography (`Nanc = pi/4u`).
``r``
    Recombination (crossover) mean rate per bp, per generation (sex-averaged), often reported in literature from direct measurement (recombination in pedigrees) or inferred from sequence data. Note you can add a crossover rate map that can modify ``r`` across the genome with ``--rec_map``. 
``u``   
    Mutation mean rate per bp, per generation, often reported in literature from mutation accumulation experiments. Note that the point mutation rate is typically used, though all mutation types with selective effects may contribute to BGS similarly, if considering different mutation types with different DFEs, see :doc:`Multiple DFEs <../guides/multiple_dfes>`.  
``g`` 
    Gene conversion initiation mean rate per bp, per generation. Note that on occasion the ``g * k`` value is reported in the literature rather than the *initiation* rate, in which case, the value should be divided by the tract length (``k``). You can add a gc initiation rate map that can modify `g` across the genome with ``--gc_map``.
``k``
    Gene conversion tract length (bp). Note that **Bvalcalc** takes only a single mean value and so does not model a distribution of tract lengths.

Distribution of fitness effects
---------------------------------

A distribution of fitness effects (DFE) describes the probability of different selective effects for new mutations when they arise. 
Selection coefficients are quantified as ``s``, the fitness effect of a mutation in homozygous state, where ``s = 0`` for a neutral allele and ``|s| = 1`` for a homozygous lethal allele. 
The effect of selection is often scaled by the effective population size (``Nanc``) as ``2*Nanc*s``. The dominance coefficient, ``h`` scales the effect of selection in heterozygous state, where ``h = 0`` for a fully recessive allele and ``h = 1`` for a fully dominant allele.

``h``
    Dominance coefficient of selected alleles. Keep at 0.5 (additive effects) unless literature suggests otherwise

When you specify a DFE in **Bvalcalc**, you are describing the probability of new deleterious mutations of a given strength arising in conserved regions (provided by e.g. the input GFF). Beneficial mutations are currently not supported.

.. note::
   **Bvalcalc** typically models a deleterious DFE consisting of non-overlapping uniform distributions between defined selection coefficient break points ranging from 0 to 1 (where 1 reflects homozygous lethal).


Basic DFE 
~~~~~~~~~~
The basic default discretized DFE contains four categories of deleterious mutations, ranging from effectively neutral (proportion described by ``f0``), to strongly deleterious (``f3``);  

See Figure 1 in `Johri et al. 2020 <https://doi.org/10.1534/genetics.119.303002>`_ for a more detailed explanation of this discretized DFE.

To specify a basic DFE, provide ``f0``, ``f1``, ``f2``, ``f3`` proportions that represent the DFE for all annotated regions in the :doc:`BED/GFF input <../introduction/bedgff_input>`. Note that the proportions must sum to 1, i.e. ``f0+f1+f2+f3 = 1``.

``f0`` 
    Proportion of effectively neutral mutations with `0 <= | 2*Nanc*s | < 1`.
    
    Note that `2*Nanc*s < 5` does not substantially contribute to BGS where HRI is not pervasive, see `Johri et al. 2020 <https://doi.org/10.1534/genetics.119.303002>`_, **Bvalcalc** will exclude the f0 proportion from BGS calculations.
``f1``
    Proportion of weakly deleterious mutations with `1 <= | 2*Nanc*s | < 10`
``f2`` 
    Proportion of moderately deleterious mutations with `10 <= | 2*Nanc*s | < 100`
``f3``
    Proportion of strongly deleterious mutations with `100 <= | 2*Nanc*s |` 

Single constant DFE 
~~~~~~~~~~~~~~~~~~~~
A fixed constant selective strength may be used for all selected sites, replacing the discretized DFE when ``--constant_dfe`` is active:

``s, proportion_synonymous``
    A single homozygous selective strength to use for all selected mutations, and the proportion of strictly neutral sites in annotated regions (e.g. synonymous).

Gamma DFE 
~~~~~~~~~~
DFE parameters may be reported in the literature as a gamma distribution. **Bvalcalc** can take gamma distribution parameters which is converted to a discretized DFE to overwrite ``f0``, ``f1``, ``f2``, ``f3`` when ``--gamma_dfe`` is active.
The gamma DFE is discretized into 9 bins with 10 break points of ``s = 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 0``. 

``mean, shape, proportion_synonymous``
    The mean and shape parameters of the gamma DFE, and the proportion of strictly neutral sites in the annotated regions (e.g. ~0.3 for synonymous sites in coding sequence annotations). 

If additional granularity is needed beyond the 9 bins, consider developing code to discretize your gamma distribution into additional bins which can be provided as input with ``--custom_dfe`` as described below.

Custom DFE 
~~~~~~~~~~

The discretized DFE can be customized with any arbitrary set of break points and proportions when ``--custom_dfe`` is active, which will overwrite the basic DFE parameters (``f0``, ``f1``, ``f2``, ``f3``). 
Note that there should be one more break point than the number of bins, the proportions should sum to 1. 

``s_breaks``
    A list of selection coefficient break points to define the bins of the discretized DFE, ranging from 0 to 1 (where 1 reflects homozygous lethal).
``bin_proportions``
    A list of proportions of new mutations in each of the bins defined as between each value in ``s_breaks``. 
    
For example, if you wanted to model a DFE with 5 bins with break points at `s = 1, 1e-2, 1e-4, 1e-6, 0` and proportions `0.1, 0.2, 0.4, 0.2, 0.1`, you would set ``s_breaks = 0, 1e-6, 1e-4, 1e-2, 1`` and ``bin_proportions = 0.1, 0.2, 0.4, 0.2, 0.1`` in the parameters file, and add the ``--custom_dfe`` flag to your CLI command.

Demography
-----------

Historical population size change as a single step-function can be accounted for with **Bvalcalc** by adding the ``--pop_change`` flag, and setting the following parameters:

``Ncur`` 
    Current population size, i.e. in the current epoch.

``time_of_change`` 
    Time in generations ago that effective population size went from ``Nanc`` to ``Ncur``. For example, if the time of change was 0.4*Nanc generations ago, and Nanc was 10000, put ``time_of_change = 4000`` or ``time_of_change = 0.4 * Nanc``.

Selfing species
---------------

Self-fertilizing (selfing) species have different evolutionary dynamics from obligate outcrossers which impact BGS (i.e., B).

The parameters file can be modified for selfing populations by altering population size (``Nanc/Ncur``), crossover rate (``r``), gene conversion rate (``g``) and dominance coefficient by an additional parameter: ``f``, which is Wright's inbreeding coefficient (F). F can be calculated from the selfing rate (F = S/(2-S)), see `Nordborg 2000 <https://doi.org/10.1093/genetics/154.2.923>`_. Note that the ``arabidopsis`` and ``pfalciparum`` default templates have the inbreeding parameter included.

For analysis of selfing populations, we recommend tailoring parameters from the ``selfing`` template. Note that the adjustments with ``f`` are coded into the template. 

.. code-block:: bash

    Bvalcalc --generate_params selfing


