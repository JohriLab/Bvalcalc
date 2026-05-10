## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## e.g. Bvalcalc --params path/to/ExampleParams.py
##
## Core parameters
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 1.33e6 / x # Ancestral population size [1]
r = 1e-8 * x # Recombination (crossover) rate per bp, per generation (sex-averaged) [2]
u = 6e-9 * x # Mutation rate (all types) per bp, per generation [1]
g = 1e-8 * x # Gene conversion initiation rate per bp, per generation [3]
k = 440 # Gene conversion tract length (bp) [3]
## Basic DFE parameters for ALL sites in annotated regions (Sum must equal 1), can overwrite with other --*_dfe flags
f0 = 0.25 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS) [4]
f1 = 0.49 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10 [4]
f2 = 0.04 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100 [4]
f3 = 0.22 # Proportion of strongly deleterious mutations with |2Ns| >= 100 [4]
## Demography parameters
Ncur = 2 * Nanc # Current population size (!Requires --pop_change) [5]
time_of_change = 0.45 * Nanc # Time in generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change) [6]
## Advanced DFE parameters 
h = 0.5 # Dominance coefficient of selected alleles [Naive value]
mean, shape, proportion_synonymous = 93.3/(2*Nanc), 0.24, 0.0 # Gamma distribution of DFE to discretize into 9 bins  [mean (s), shape, strictly neutral proportion] (!Requires --gamma_dfe) [7]
s_breaks = 0, 1/(2*Nanc), 10/(2*Nanc), 100/(2*Nanc), 1 # Custom DFE parameter controlling the homozygous selection coefficient (s) breakpoints (!Requires --custom_dfe) [Naive value]
bin_proportions = 0.25, 0.25, 0.25, 0.25 # Custom DFE parameter controlling the proportion of mutations between each bin by s_breaks, overwriting the f1-f3 values above (!Requires --custom_dfe) [Naive value]
## Literature cited
# [1] Keightley et al 2014  doi: 10.1534/genetics.113.158758
# [2] Comeron et al 2012 doi: 10.1371/journal.pgen.1002905
# [3] Miller et al 2016 doi: 10.1534/genetics.115.186486
# [4] Johri et al 2020 doi: 10.1534/genetics.119.303002
# [5] Laurent et al 2011 doi: 10.1093/molbev/msr031
# [6] Kapopoulou et al 2018 doi: 10.1093/gbe/evy185
# [7] Daigle et al. 2026 doi: 10.64898/2026.03.01.708907v1