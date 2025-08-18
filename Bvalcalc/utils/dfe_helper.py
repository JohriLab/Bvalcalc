import scipy.stats as st
import os
import importlib.util
from typing import Dict, Any
from functools import lru_cache

GAMMA_DFE = False  # Default, instead of prop injected
CONSTANT_DFE = False # Default, instead of prop injected

def get_DFE_params(params_path: str | None = None, gamma_dfe: bool = False, constant_dfe: bool = False) -> Dict[str, Any]:
    """
    Load and validate population parameters from the file pointed to by
    `params_path` or, if None, by the BCALC_params env var.
    Returns a dictionary of parameters for use in B-value calculations.
    """
    # 1. Determine the path: either passed in or from the env var
    if params_path is None:
        params_path = os.environ.get("BCALC_params")
        if not params_path:
            raise KeyError(
                "Environment variable BCALC_params not set. "
                "Cannot load pop-gen parameters."
            )

    # 2. Load the module
    spec = importlib.util.spec_from_file_location("params", params_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load spec for params from {params_path}")
    pop = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pop)

    # 3. Extract and validate required attributes
    required_names = ["g", "k", "r", "u", "Nanc", "h", "f0", "f1", "f2", "f3"]
    params: Dict[str, Any] = {}
    for name in required_names:
        if not hasattr(pop, name):
            raise AttributeError(f"params must define '{name}'")
        val = getattr(pop, name)
        if val is None:
            raise AttributeError(f"params parameter '{name}' is None; must be a numeric value")
        params[name] = float(val)

    # 4. Optional gamma-DFE override
    if GAMMA_DFE or gamma_dfe is not False: # The GAMMA_DFE is prop injected by CLI, gamma_dfe is provided by API
        mean = getattr(pop, 'mean', None)
        shape = getattr(pop, 'shape', None)
        prop_syn = getattr(pop, 'proportion_synonymous', None)
        if mean is None or shape is None or prop_syn is None:
            raise AttributeError(
                "params must define 'mean', 'shape' and 'proportion_synonymous' when --gamma_dfe is active"
            )
        from .dfe_helper import gammaDFE_to_discretized
        f0, f1, f2, f3 = gammaDFE_to_discretized(mean, shape, prop_syn)
        params.update({"f0": f0, "f1": f1, "f2": f2, "f3": f3})


    # 5. Set derived parameters and thresholds
    params["gamma_cutoff"] = 5
    params["t0"] = 0.0
    Nanc = params["Nanc"]
    h = params["h"]
    # Calculate generation-scale thresholds
    params["t1"] = h * (1.0 / (2.0 * Nanc))
    params["t1half"] = h * (params["gamma_cutoff"] / (2.0 * Nanc))
    params["t2"] = h * (10.0 / (2.0 * Nanc))
    params["t3"] = h * (100.0 / (2.0 * Nanc))
    params["t4"] = h * 1.0

    if CONSTANT_DFE or constant_dfe is not False: # The CONSTANT_DFE is prop injected by CLI, constant_dfe is provided by API
        s = getattr(pop, "s", None)
        prop_syn = getattr(pop, 'proportion_synonymous', None)
        if s is None or prop_syn is None:
            raise AttributeError(
                "params must define 's' and 'proportion_synonymous' when --constant_dfe is active"
            )
        params["t_constant"] = h * s # Set parameter to be exported to calculateB
    else:         
        params["t_constant"] = None

    return params


def gammaDFE_to_discretized(mean: float, shape: float, proportion_synonymous: float):
    if mean <= 0 or shape <= 0:
        raise ValueError("`mean` and `shape` must be positive.")
    if not (0 <= proportion_synonymous < 1):
        raise ValueError("`proportion_synonymous` must be in [0, 1).")

    theta = mean / shape               # scale parameter
    dist  = st.gamma(a=shape, scale=theta)

    # cumulative probabilities at the cut‑points
    c1   = dist.cdf(1.0)
    c10  = dist.cdf(10.0)
    c100 = dist.cdf(100.0)


    f0 = c1                      # (0, 1]
    f1 = c10  - c1               # (1,10]
    f2 = c100 - c10              # (10,100]
    f3 = 1.0   - c100            # >100

    # 4) scale to sum to (1 - p_syn)
    scale = 1.0 - proportion_synonymous
    f0, f1, f2, f3 = (f0 * scale,
                      f1 * scale,
                      f2 * scale,
                      f3 * scale)

    # 5) add synonymous fraction back into f0
    f0 += proportion_synonymous

    print(f"Converting gamma distribution to discretized DFE")
    print(f"Gamma params: mean = {mean}, shape = {shape}, scale = {theta}")
    print(f"Inferred f0, f1, f2, f3 = ", f0, f1, f2, f3)

    return f0, f1, f2, f3