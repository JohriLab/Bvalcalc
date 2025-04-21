## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## Usage: ./Bvalcalc --region --pop_params ExampleParams.py
##        ./Bvalcalc --genome --pop_params ExampleParams.py
##
## Core parameters
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 1e6/x # Ancestral population size
r = 1*1e-8*x # Recombination (crossover) rate per bp, per generation
u = 3*1e-9*x # Mutation rate (all types) per bp, per generation
g = 1*1e-8*x # Gene conversion initiation rate per bp, per generation
k = 440 # Gene conversion tract length (bp)
##
## Demography parameters
Ncur = Nanc*2 # Current population size (!Requires --pop_change)
time_of_change = 0.1 # Time in 2Ncur generations ago that effective population size went from Nanc to Ncur (!Requires --pop_change)
##
## Distribution of fitness effects (DFE) parameters
f0 = 0.25 # Proportion of effectively neutral mutations with 0 <= |2Nes| < 1 (e.g.)
f1 = 0.49 # Proportion of weakly deleterious mutations with 1 <= |2Nes| < 10
f2 = 0.04 # Proportion of moderately deleterious mutations with 10 <= |2Nes| < 100
f3 = 0.22 # Proportion of strongly deleterious mutations with |2Nes| >= 100
##
## Advanced DFE parameters 
h=0.5 # Dominance coefficient of selected alleles
##
mean, shape = 500, 0.5 # Gamma distribution of DFE [mean in 2Ns, shape parameter] (!Requires --gamma_dfe)
##
## !DO NOT CHANGE BELOW UNLESS WELL INFORMED BY POPGEN THEORY!
gamma_cutoff = 5 # 2Ns threshold for effectively neutral alleles, mutations below this threshold will be ignored in B calculation. Keep as 5 unless theory suggests otherwise.
t0 = 0.0 # Start of neutral class (t=hs=0)
t1 = h*(1/(2*Nanc)) # Start of f1 class (2Ns=1)
t1half = h*(gamma_cutoff/(2*Nanc)) # 2Ns threshold for effectively neutral alleles
t2 = h*(10/(2*Nanc)) # End of f1 class, start of f2 class (2Ns=10)
t3 = h*(100/(2*Nanc)) # End of f2 class, start of f3 class (2Ns=100) 
t4 = h*1.0 # End of f3 class (s=1)