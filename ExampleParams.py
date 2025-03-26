## Population genetic parameters for the simulated or empirical population
## Accurate estimation requires accurate and appropriate parameters
##
## Usage: ./Bvalcalc --region --pop_params ExampleParams.py
##        ./Bvalcalc --genome --pop_params ExampleParams.py

x = 1 # Scaling Factor (N,u,r), keep as [1] unless calculating for rescaled simulations
g = 1*1e-8*x # Gene conversion initiation rate per bp, per generation
k = 440 # Gene conversion tract length (bp)
r = 1*1e-8*x #rate of recombination
u = 3*1e-9*x #(*Mutation rate*)

Nanc = 1e6/x #(Ancestral population size)
Ncur = 1e6/x #(Current population size)
gamma_cutoff = 5 #5.0
h=0.5
t0 = 0.0
t1 = h*(1/(2*Nanc))
t1half = h*(gamma_cutoff/(2*Nanc)) #(* This is the cut-off value of 2Nes=5. This derivation assumes that all mutations with 2Nes<5 will not contribute to BGS *)
t2 = h*(10/(2*Nanc))
t3 = h*(100/(2*Nanc))
t4 = h*1.0

f0 = 0.25 #(*Proportion of effectively neutral mutations with 0 <= |2Nes| < 1 *) 0.25
f1 = 0.49 #(*Proportion of weakly deleterious mutations with 1 <= |2Nes| < 10 *) 0.6533 0.49
f2 = 0.04 #(*Proportion of moderately deleterious mutations with 10 <= |2Nes| < 100 *) 0.0533 0.04
f3 = 0.22 #(*Proportion of strongly deleterious mutations with |2Nes| >= 100 *) 0.2933 0.22

time_of_change = 1 #0.1/0.5/1(This is the time of change in 2Ncur  in the past generations.)