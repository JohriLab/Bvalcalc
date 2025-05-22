Bvalcalc Documentation
======================

Welcome to the documentation for **Bvalcalc**, a CLI tool for easily calculating B-values and B-maps. 

Background selection (BGS) is the process by which diversity is reduced at neutral sites linked to conserved regions such as exons that experience direct purifying selection, and across the genome. **Bvalcalc** allows users to easily calculate the relative reduction of diversity due to BGS (B) using several different modes. Calculating B-values is important for understanding the linked effects of selection, can be used to establish an evolutionary null model for inference, and can aid in selecting the most neutrally evolving sites, see :doc:`Demographic Inference with B-map <guides/demography>`. 

To learn more about background selection, see :doc:`About BGS and B-values <prebuilt_maps/learn>`. To access pre-built B-maps for select species, see :doc:`Query Species B-maps <prebuilt_maps/query>`. To calculate B-values for yourself, install **Bvalcalc** and explore the docs, we recommend getting started with the :doc:`Quickstart Tutorial <guides/quickstart>`.

.. code-block:: console

    Bvalcalc [--generate_params [SPECIES]|--site|--gene|--region|--genome]

   --generate_params [SPECIES]
                         Save population parameters from a species template (human, drosophila, arabidopsis, celegans, mouse, pfalciparum, selfing)
   --site, -s            Calculate B values for a single site from a selected element
   --gene, -g            Calculate B values for a region adjacent to a single selected element
   --region, -r          Calculate B values for a chromosomal region, considering genome-wide effects
   --genome, -w          Calculate B values genome-wide for all sites considering all elements

.. toctree::
   :maxdepth: 1
   :caption: Introduction

   introduction/installation
   introduction/bedgff_input
   introduction/generate_params

.. toctree::
   :maxdepth: 1
   :caption: Modules

   modes/site
   modes/gene
   modes/region
   modes/genome
   modes/vcf

.. toctree::
   :maxdepth: 2
   :caption: Guides

   guides/quickstart
   guides/params
   guides/multiple_dfes
   guides/demography
   guides/sweepfinder2

.. toctree::
   :maxdepth: 1
   :caption: API

   api/calculateB

.. toctree::
   :hidden:
   :caption: About BGS and B-values
   :maxdepth: 1

   prebuilt_maps/learn

.. toctree::
   :hidden:
   :caption: Query species B-maps
   :maxdepth: 1

   prebuilt_maps/query