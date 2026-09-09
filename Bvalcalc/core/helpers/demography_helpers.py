import os
import numpy as np
import importlib.util

def _load_params_module():
    """
    Dynamically load the population parameters from BCALC_params.
    """
    params_path = os.environ.get("BCALC_params")
    if not params_path:
        raise KeyError(
            "Environment variable BCALC_params not set. Cannot load population parameters."
        )
    spec = importlib.util.spec_from_file_location("params", params_path)
    params = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(params)
    return params

def get_Bcur(Banc: np.ndarray) -> np.ndarray:
    """
    Apply demographic change to an array of B values using parameters from the population parameters module.

    Args:
        Banc: np.ndarray of B under ancestral population size

    Returns:
        np.ndarray of B under current population size
    """
    pop = _load_params_module()
    Nanc = pop.Nanc
    Ncur = pop.Ncur
    time_of_change = pop.time_of_change

    # Ratio of ancestral to current population sizes
    R = Nanc / Ncur

    # Compute exponential terms
    exp_term_num = np.exp(-(time_of_change / (2 * Ncur)) / Banc)
    exp_term_den = np.exp(-(time_of_change / (2 * Ncur)))

    # Numerator and denominator of the demographic adjustment formula
    numerator = Banc * (1 + (R - 1) * exp_term_num)
    denominator = 1 + (R - 1) * exp_term_den

    if Ncur < Nanc:
        print(f"WARNING: Calculating BGS under a population contraction with factor Ncur/Nanc = {1/R}. This will typically increase genetic drift and lead to weaker BGS effects for mutations with selection coefficients that fall below 2*Ncur*s ~ 5. The --pop_change option does not currently account for the shift in these mutation's BGS effects toward more neutral dynamics, so calculated B-values may be lower than expected near conserved regions. Check the proportion of mutations that will transition from -2Ns > 5 to -2Ns < 5 in the contracted population to assess the severity of BGS overestimation.")

    if Nanc > Ncur:
        print(f"WARNING: Calculating BGS under a population expansion with factor Ncur/Nanc = {1/R}. For mutations with selection coefficients with -2*Nanc*s < 5 but -2*Ncur*s > 5, the --pop_change option does not currently account for the shift in these mutation's BGS effects toward more strongly selected dynamics, so calculated B-values may be higher than expected near conserved regions. This is typically a small proportion of the DFE and is unlikely to have strong effects on most analyses unless a large proportion of selected mutations transition from -2Ns < ~5 to -2Ns > ~5.")

    return numerator / denominator