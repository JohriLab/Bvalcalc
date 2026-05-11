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
## Basic DFE parameters for ALL sites in annotated regions (Sum must equal 1), can overwrite with other --*_dfe flags
f0 = 0.25 # Proportion of effectively neutral mutations with 0 <= |2Ns| < 1 (Note that 2Ns<5 does not contribute to BGS) [Naive value]
f1 = 0.25 # Proportion of weakly deleterious mutations with 1 <= |2Ns| < 10 [Naive value]
f2 = 0.25 # Proportion of moderately deleterious mutations with 10 <= |2Ns| < 100 [Naive value]
f3 = 0.25 # Proportion of strongly deleterious mutations with |2Ns| >= 100 [Naive value]
# ## Demography parameters
Ncur = Nanc # Current population size (!Requires --pop_change) [Naive value] 
time_of_change = 1 * Nanc # Time in generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change) [Naive value]
## Advanced DFE parameters 
h = 0.5 + (f-0.5*f) # Dominance coefficient of selected alleles, NOTE: this is h_eff for h=0.5, replace BOTH 0.5's with your dominance coefficient [Naive value]
mean, shape, proportion_synonymous = 100 / (2*Nanc), 1, 0.3 # Gamma distribution of DFE to discretize into 9 bins [mean (s), shape, strictly neutral proportion] (!Requires --gamma_dfe) [Naive value]
s_breaks = 0, 1/(2*Nanc), 10/(2*Nanc), 100/(2*Nanc), 1 # Custom DFE parameter controlling the homozygous selection coefficient (s) breakpoints (!Requires --custom_dfe) [Naive value]
bin_proportions = 0.25, 0.25, 0.25, 0.25 # Custom DFE parameter controlling the proportion of mutations between each bin by s_breaks, overwriting the f1-f3 values above (!Requires --custom_dfe) [Naive value]
## Literature cited
# [1]
# [2]
# [3]
# [4]
