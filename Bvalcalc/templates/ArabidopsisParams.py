## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## e.g. Bvalcalc --params path/to/ExampleParams.py
##
## Core parameters
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 3.9e5 / x # Ancestral population size [d]
r = 7.465e-8 * x # Recombination (crossover) rate per bp, per generation (sex-averaged) [a]
u = 6.95e-9 * x # Mutation rate (all types) per bp, per generation [b]
g = r / 140 * x # Gene conversion initiation rate per bp, per generation [c]
k = 50 # Gene conversion tract length (bp) [c]
## DFE parameters for ALL sites in annotated regions (Sum must equal 1)
Xf0 = 0.25 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS) [4]
Xf1 = 0.49 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10 [4]
Xf2 = 0.04 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100 [4]
Xf3 = 0.22 # Proportion of strongly deleterious mutations with |2Ns| >= 100 [4]
## Demography parameters
XNcur = 2 * Nanc # Current population size (!Requires --pop_change) [5]
Xtime_of_change = 0.45 # Time in Nanc generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change) [6]
## Advanced DFE parameters 
Xh = 0.5 # Dominance coefficient of selected alleles [Naive value]
Xmean, shape, proportion_synonymous = 500, 0.5, 0.3 # Gamma distribution of DFE to discretize and replace f0-f3 [mean (2Ns), shape, proportion synonymous] (!Requires --gamma_dfe) [Naive value]## Literature cited
# [a] Rowan et al. 2019 doi: 10.1534/genetics.119.302406
# [b] Weng et al. 2018 doi: 10.1534/genetics.118.301721
# [c] Wijnker et al. 2013 doi: 10.7554/eLife.01426
# [d] Nordborg et al. 2005 doi: 10.1371/journal.pbio.0030196