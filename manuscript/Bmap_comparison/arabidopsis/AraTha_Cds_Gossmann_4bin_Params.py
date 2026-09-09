## Arabidopsis thaliana CDS parameters with the Gossmann DFE in four 2*Nanc*s classes

## Core parameters retained from the exact colleague-template run
f = 0.97
x = 1
Nanc = 1e5 / (1 + f) / x
r = 7.465e-8 * (1 - f) * x
u = 6.95e-9 * x
g = r * 50 * (1 - f) * x
k = 553

## Gamma(mean 2*Nanc*s=3950, shape=0.177) integrated over the four
## legacy 2*Nanc*s classes. The 30% synonymous mass is included in f0.
f0 = 0.428656619
f1 = 0.064720817
f2 = 0.097120240
f3 = 0.409502324

## Demography retained from the exact colleague-template run
Ncur = 0.5 * Nanc
time_of_change = 1 * Nanc

h = 0.5 + (f - 0.5 * f)

