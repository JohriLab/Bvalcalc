import os
import numpy as np
import importlib.util

def _load_pop_params_module():
    """
    Dynamically load the population parameters from BCALC_POP_PARAMS.
    """
    params_path = os.environ.get("BCALC_POP_PARAMS")
    if not params_path:
        raise KeyError(
            "Environment variable BCALC_POP_PARAMS not set. Cannot load population parameters."
        )
    spec = importlib.util.spec_from_file_location("pop_params", params_path)
    pop_params = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pop_params)
    return pop_params

def get_Bcur(Banc: np.ndarray) -> np.ndarray:
    """
    Apply demographic change to an array of B values using parameters from the population parameters module.

    Args:
        Banc: np.ndarray of B under ancestral population size

    Returns:
        np.ndarray of B under current population size
    """
    pop = _load_pop_params_module()
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

    return numerator / denominator
