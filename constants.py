# Constants (make sure to define these in your script or pass them as needed)
x = 100 # Scaling factor to rescale N, u, r
g = 1*1e-8*x #1e-8 #rate of gene conversion
tract_len=440 #mean tract length of gene conversion in base pairs
r = 0.5*1e-8*x #rate of recombination
u = 3*1e-9*x #(*Mutation rate*)

Nanc = 1e6/x
Ncur = 1e6/x #(Current population size)
gamma_cutoff = 5 #5.0
h=0.5
t0 = 0.0
t1 = h*(1/(2*Nanc))
t1half = h*(gamma_cutoff/(2*Nanc)) #(* This is the cut-off value of 2Nes=5. This derivation assumes that all mutations with 2Nes<5 will not contribute to BGS *)
t2 = h*(10/(2*Nanc))
t3 = h*(100/(2*Nanc))
t4 = h*1.0

f0 = 0.1 #(*Proportion of effectively neutral mutations with 0 <= |2Nes| < 1 *) 0.25
f1 = 0.2 #(*Proportion of weakly deleterious mutations with 1 <= |2Nes| < 10 *) 0.6533 0.49
f2 = 0.3 #(*Proportion of moderately deleterious mutations with 10 <= |2Nes| < 100 *) 0.0533 0.04
f3 = 0.4 #(*Proportion of strongly deleterious mutations with |2Nes| >= 100 *) 0.2933 0.22

time_of_change=1.0 #0.1/0.5/1(This is the time of change in 2Ncur  in the past generations.)