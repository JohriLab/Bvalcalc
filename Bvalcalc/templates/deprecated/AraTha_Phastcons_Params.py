## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## e.g. Bvalcalc --params path/to/ExampleParams.py
##
## Core parameters
f = 0.97 # Selfing rate (f = S/(2-S); Wright's inbreeding coefficient) [1]
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 1e5 / (1+f) / x # Ancestral population size [2]
r = 7.465e-8 * (1-f) * x # Recombination (crossover) rate per bp, per generation (sex-averaged) [3]
u = 6.95e-9 * x # Mutation rate (all types) per bp, per generation [4]
g = r * 50 * (1-f) * x # Gene conversion initiation rate per bp, per generation [5]
k = 553 # Gene conversion tract length (bp) [5]
## DFE parameters for ALL sites in annotated regions (Sum must equal 1)
f0 = 0.28 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS) [6*]
f1 = 0.33 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10 [6*]
f2 = 0.35 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100 [6*]
f3 = 0.04 # Proportion of strongly deleterious mutations with |2Ns| >= 100 [6*]
## Demography parameters
Ncur = 0.5 * Nanc # Current population size (!Requires --pop_change) [2]
time_of_change = 1 # Time in Nanc generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change) [2]
## Advanced DFE parameters 
h = 0.5 + (f-0.5*f) # Dominance coefficient of selected alleles [Naive value]
mean, shape, proportion_synonymous = 500, 0.5, 0.3 # Gamma distribution of DFE to discretize and replace f0-f3 [mean (2Ns), shape, proportion synonymous] (!Requires --gamma_dfe) [Naive value]## Literature cited
# [1] Platt et al. 2010 doi: 10.1371/journal.pgen.1000843
# [2] Durvasula et al. 2017 doi: 10.1073/pnas.1616736114
# [3] Rowan et al. 2019 doi: 10.1534/genetics.119.302406
# [4] Weng et al. 2018 doi: 10.1534/genetics.118.301721
# [5] Yang et al. 2012 doi: 10.1073/pnas.1211827110
# [6] Williamson et al 2014 doi: 10.1371/journal.pgen.1004622 *Note that this is an estimate for phastCons results for Capsella rubella as A. thaliana results are not available