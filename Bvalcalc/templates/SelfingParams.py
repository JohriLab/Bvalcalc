## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## e.g. Bvalcalc --params path/to/TemplateParams.py
##
## Core parameters
f = 0 # Selfing rate (f = S/(2-S); Wright's inbreeding coefficient) [Naive value]
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 1e5 / (1+f) / x # Ancestral population size [Naive value]
r = 1e-8 * (1-f) * x # Recombination (crossover) rate per bp, per generation (sex-averaged) [Naive value]
u = 1e-8 * x # Mutation rate (all types) per bp, per generation [Naive value]
g = 1e-8 * (1-f) * x # Gene conversion initiation rate per bp, per generation [Naive value]
k = 500 # Gene conversion tract length (bp) [Naive value]
## DFE parameters for ALL sites in annotated regions (Sum must equal 1)
f0 = 0.25 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS) [Naive value]
f1 = 0.25 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10 [Naive value]
f2 = 0.25 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100 [Naive value]
f3 = 0.25 # Proportion of strongly deleterious mutations with |2Ns| >= 100 [Naive value]
# ## Demography parameters
Ncur = Nanc # Current population size (!Requires --pop_change) [Naive value] 
time_of_change = 1 # Time in Nanc generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change) [Naive value]
## Advanced DFE parameters 
h = 0.5 + (f-0.5*f) # Dominance coefficient of selected alleles [Naive value]
mean, shape, proportion_synonymous = 100, 1, 0.3 # Gamma distribution of DFE to discretize and replace f0-f3 [mean (2Ns), shape, proportion synonymous] (!Requires --gamma_dfe) [Naive value]
## Literature cited
# [1]
# [2]
# [3]
# [4]
