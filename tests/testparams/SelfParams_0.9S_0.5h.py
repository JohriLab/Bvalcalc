## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## Usage: ./Bvalcalc --region --pop_params ExampleParams.py
##        ./Bvalcalc --genome --pop_params ExampleParams.py
## poetry run Bvalcalc --gene --pop_params tests/testparams/SelfParams_0.9S_0.5h.py

## Core parameters
f = 0.9 # Selfing rate (F = S/(2-S); Wright's inbreeding coefficient)
# S = 0.9473684211
x = 100 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 1e6 / (1+f) /x # Ancestral population size
r = 0.5*1e-8* (1-f) *x # Recombination (crossover) rate per bp, per generation
u = 3*1e-9*x # Mutation rate (all types) per bp, per generation
g = 0*1e-8 * (1-f) *x # Gene conversion initiation rate per bp, per generation
k = 440 # Gene conversion tract length (bp)

## Demography parameters
Ncur = Nanc # Current population size (!Requires --pop_change)
time_of_change = 1 # Time in 2Ncur generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change)

# Distribution of fitness effects (DFE) parameters (Must equal 1)
f0 = 0.1 #(*Proportion of effectively neutral mutations with 0 <= |2Nes| < 1 *) 0.25
f1 = 0.2 #(*Proportion of weakly deleterious mutations with 1 <= |2Nes| < 10 *) 0.6533 0.49
f2 = 0.3 #(*Proportion of moderately deleterious mutations with 10 <= |2Nes| < 100 *) 0.0533 0.04
f3 = 0.4 #(*Proportion of strongly deleterious mutations with |2Nes| >= 100 *) 0.2933 0.22

## Advanced DFE parameters
gamma_cutoff = 5 # 2Ns threshold for effectively neutral alleles, mutations below this threshold will be ignored in B calculation. Keep as 5 unless theory suggests otherwise.
h=0.5 + (f-0.5*f) # Dominance coefficient of selected alleles
mean, shape, proportion_synonymous = 500, 0.5, 0.3 # Gamma distribution of DFE to discretize and replace f0-f3 [mean (2Ns), shape] (!Requires --gamma_dfe)