"""D. melanogaster CDS parameters used for the paper-map v1.4.4 rerun."""

x = 1
Nanc = 1.225393e6 / x
r = 1e-8 * x
u = 3e-9 * x
g = 1e-8 * x
k = 440

# Required legacy values. These are replaced when --gamma_dfe is active.
f0, f1, f2, f3 = 0.25, 0.49, 0.04, 0.22

Ncur = 1.10802003929 * Nanc
time_of_change = 1
h = 0.5

mean = 811 / (2 * Nanc)
shape = 0.347
proportion_synonymous = 0.3

