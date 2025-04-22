## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## Usage: ./Bvalcalc --region --pop_params ExampleParams.py
##        ./Bvalcalc --genome --pop_params ExampleParams.py
##
## Core parameters
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 1e6/x # Ancestral population size
r = 1e-8*x # Recombination (crossover) rate per bp, per generation
u = 3e-9*x # Mutation rate (all types) per bp, per generation
g = 1e-8*x # Gene conversion initiation rate per bp, per generation
k = 440 # Gene conversion tract length (bp)
## Demography parameters
Ncur = Nanc*2 # Current population size (!Requires --pop_change)
time_of_change = 0.1 # Time in 2Ncur generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change)
## Distribution of fitness effects (DFE) parameters
f0 = 0.25 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS)
f1 = 0.49 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10
f2 = 0.04 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100
f3 = 0.22 # Proportion of strongly deleterious mutations with |2Ns| >= 100
##
## Advanced DFE parameters 
h=0.5 # Dominance coefficient of selected alleles
mean, shape = 500, 0.5 # Gamma distribution of DFE to replace f0-f3 [mean (2Ns), shape] (!Requires --gamma_dfe)
##

## !DO NOT CHANGE BELOW UNLESS WELL INFORMED BY POPGEN THEORY!
