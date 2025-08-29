## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## e.g. Bvalcalc --params path/to/ExampleParams.py
##
## Core parameters
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 7300 / x # Ancestral population size [1]
r = 1e-8 * x # Recombination (crossover) rate per bp, per generation (sex-averaged) [2]
u = 1.25e-8 * x # Mutation rate (all types) per bp, per generation [2]
g = 3.27e-8 * x # Gene conversion initiation rate per bp, per generation [3]
k = 113 # Gene conversion tract length (bp) [3]
## DFE parameters for ALL sites in annotated regions (Sum must equal 1)
f0 = 0.51 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS) [4]
f1 = 0.14 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10 [4]
f2 = 0.14 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100 [4]
f3 = 0.21 # Proportion of strongly deleterious mutations with |2Ns| >= 100 [4]
## Demography parameters
Ncur = 14474 # Current population size (!Requires --pop_change) [1]
time_of_change = 0.81 # Time in Nanc generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change) [1]
## Advanced DFE parameters 
h = 0.5 # Dominance coefficient of selected alleles [Naive value]
mean, shape = 500, 0.5 # Gamma distribution of DFE to discretize and replace f0-f3 [mean (2Ns), shape] (!Requires --gamma_dfe) [Naive value]
## Literature cited
# [1] Gutenkunst et al 2009 doi: 10.1371/journal.pgen.1000695
# [2] Kong et al 2012 doi: 10.1038/nature11396
# [3] Palsson 2025 doi: 10.1038/s41586-024-08450-5
# [4] Huber et al 2017 doi: 10.1073/pnas.1619508114