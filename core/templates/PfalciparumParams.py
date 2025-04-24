## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## e.g. Bvalcalc --pop_params path/to/ExampleParams.py
##
## Core parameters
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 5e4 / x # Ancestral population size (0.9)
r = 7.4e-7 * x # Recombination (crossover) rate per bp, per generation (sex-averaged)
u = 5.6e-9 * x # Mutation rate (all types) per bp, per generation
# ?g = 1e-8 * x # Gene conversion initiation rate per bp, per generation
# ?k = 440 # Gene conversion tract length (bp)
# ## DFE parameters (Sum must equal 1)
f0 = 0.4 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS) [Naive value]
f1 = 0.2 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10 [Naive value]
f2 = 0.2 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100 [Naive value]
f3 = 0.2 # Proportion of strongly deleterious mutations with |2Ns| >= 100 [Naive value]
# ## Demography parameters
# ?Ncur = Nanc*2 # Current population size (!Requires --pop_change)
# ?time_of_change = 0.1 # Time in Nanc generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change)
## Advanced DFE parameters 
h = 0.5 # Dominance coefficient of selected alleles
mean, shape = 500, 0.5 # Gamma distribution of DFE to discretize and replace f0-f3 [mean (2Ns), shape] (!Requires --gamma_dfe)
## Literature cited
# [1]
# [2]
# [3]
# [4]
