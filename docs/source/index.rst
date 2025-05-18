Bvalcalc Documentation
======================

Welcome to the documentation for **Bvalcalc**, a CLI tool for easily calculating B-values and B-maps. 

Background selection (BGS) is the process by which diversity is reduced at neutral sites linked to conserved regions such as exons that experience direct purifying selection, and across the genome. **Bvalcalc** allows users to easily calculate the relative reduction of diversity due to BGS (B) using several different modes. Calculating B-values is important for understanding the linked effects of selection, can be used to establish an evolutionary null model for inference, and can aid in selecting the most neutrally evolving sites for demographic inference. 

To learn more about background selection, see [About BGS and B-values]. To access pre-built B-maps for select species, see [Query Species B-maps]. To calculate B-values for yourself, install **Bvalcalc** and explore the docs, we recommend getting started with the [Quickstart Tutorial].


.. toctree::
   :maxdepth: 1
   :caption: Introduction

   introduction/cli_overview
   introduction/installation
   introduction/generate_params

.. toctree::
   :maxdepth: 1
   :caption: Modules

   modes/site
   modes/gene
   modes/region
   modes/genome

.. toctree::
   :maxdepth: 1
   :caption: Guides

   guides/quickstart
   guides/params
   guides/selfing
   guides/demography
   guides/multiple_dfes

.. toctree::
   :maxdepth: 1
   :caption: API

   api

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