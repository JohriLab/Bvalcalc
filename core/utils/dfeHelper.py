import sys
import scipy.stats as st
import os, importlib.util
from typing import Tuple

def getDFEparams(gamma_dfe = None) -> Tuple[
        float, float, float, float,   # g, k, r, u
        float, float,                 # Nanc, h
        float, float, float, float,   # f0‑f3
        float, float,                 # gamma_cutoff, t0
        float, float, float, float, float]:   # t1, t1half, t2, t3, t4
    """
    Load pop‑gen parameters from the params file pointed to by $BCALC_POP_PARAMS
    and return everything needed by calculateB.py as a single tuple.
    """

    spec = importlib.util.spec_from_file_location("pop_params",
                                                  os.environ["BCALC_POP_PARAMS"])
    _pop = importlib.util.module_from_spec(spec); spec.loader.exec_module(_pop)
    g, k, r, u, Nanc, h, f0, f1, f2, f3, mean, shape = (
        getattr(_pop, v) for v in
        ['g', 'k', 'r', 'u', 'Nanc', 'h', 'f0', 'f1', 'f2', 'f3', 'mean', 'shape']
    )

    if gamma_dfe is not None:
        f0, f1, f2, f3 = gammaDFE_to_discretized(mean, shape)

    gamma_cutoff = 5 # 2Ns threshold for effectively neutral alleles, mutations below this threshold will be ignored in B calculation. Keep as 5 unless theory suggests otherwise.
    t0 = 0.0 # Start of neutral class (t=hs=0)
    t1 = h*(1/(2*Nanc)) # Start of f1 class (2Ns=1)
    t1half = h*(gamma_cutoff/(2*Nanc)) # 2Ns threshold for effectively neutral alleles
    t2 = h*(10/(2*Nanc)) # End of f1 class, start of f2 class (2Ns=10)
    t3 = h*(100/(2*Nanc)) # End of f2 class, start of f3 class (2Ns=100) 
    t4 = h*1.0 # End of f3 class (s=1)

    return (g, k, r, u, Nanc, h, f0, f1, f2, f3,
            gamma_cutoff, t0, t1, t1half, t2, t3, t4)


def gammaDFE_to_discretized(mean: float, shape: float):
    if mean <= 0 or shape <= 0:
        raise ValueError("Both mean and shape must be positive.")

    theta = mean / shape               # scale parameter
    dist  = st.gamma(a=shape, scale=theta)

    # cumulative probabilities at the cut‑points
    c1   = dist.cdf(1.0)
    c10  = dist.cdf(10.0)
    c100 = dist.cdf(100.0)

    f0 = c1                      # (0,1]
    f1 = c10  - c1               # (1,10]
    f2 = c100 - c10              # (10,100]
    f3 = 1.0   - c100            # >100

    print(f"Converting gamma distribution to discretized DFE")
    print(f"Gamma params: mean = {mean}, shape = {shape}, scale = {theta}")
    print(f"Inferred f0, f1, f2, f3 = ", f0, f1, f2, f3)

    return f0, f1, f2, f3
