# Bvalcalc CLI Quickstart

## Installation

You can install Bvalcalc via pip:

```bash
pip install bvalcalc
```

## Specifying popgen parameters

First we need to copy population genetic parameters from one of the templates, e.g. Drosophila.
In your own analysis you'll need to check the literature and use informed values for your population.

```bash
./Bvalcalc.py --generate_params drosophila
```

## Calculating B recovery from a single element

Now, let's calculate the recovery of B as a function of distance from a single conserved element [10kb by default] and save it to gene_B.png

```bash
./Bvalcalc.py \
--gene \
--pop_params ./DrosophilaParams.py \
--plot_output gene_B.png
```

## Calculating B for a region of the genome

Next, we can use a BED file that specifies which genomic ranges are conserved to calculate B for a region in the genome.
Let's calculate B for a 1Mb region in the middle of chromosome 2R [9500000-10500000] and save it to 1Mb_B.png.

```bash
./Bvalcalc.py \
--region \
--calc_region 2R:9500000-10500000 \
--pop_params ./DrosophilaParams.py \
--bedgff_path examples/dmel6_2R_genes.csv \
--plot_output 1Mb_B.png
```

Have a look at the plots, the blue sections of the graph indicate neutral regions and black indicates conserved elements.
That's all that's necessary for many analyses, especially if you're only interested in B values for a specific region of the genome, or are testing against simulated results.

## Calculating a complete B-map

If you wanted to generate a complete B-map for all sites across all chromosomes you can use the following command, though note its a lot more data to crunch!
The following is the basic command, recommended options include `--chr_sizes` with a file specifying each chromosome's length, `--pop_change` for demography, `--rec_map` and `--gc_map` for crossover and gene conversion maps.

```bash
./Bvalcalc.py \
--genome \
--pop_params ./DrosophilaParams.py \
--bedgff_path examples/dmel6_2R_genes.csv
```

A caveat to the `--region` and `--genome` modes is that by default they combine and simplify distant elements in discrete chunks which can slightly change the distance of distant conserved elements when
calculating B. The default chunk size is 20kb and the window within which calculations are perfectly precise is three chunks in each direction (140kb total). This allows for vastly improved performance
and typically will not result in directional biases of B estimates. However, to achieve more exact results you can specify the size of the chunks with `--chunk_size`, and the size of the window to
perform perfectly precise calculations with `--precise_chunks`, though this will come at the cost of perfomance so consider using HPC resources or limiting to a specific region with `--region`.
