import importlib.util
import numpy as np
import os

# Internal cache
_pop_params = None

def get_pop_params():
    global _pop_params
    if _pop_params is None:
        # Read the path from an environment variable
        path = os.environ.get("BCALC_POP_PARAMS")
        if path is None:
            raise RuntimeError("Environment variable BCALC_POP_PARAMS not set — please provide a path to your parameter file.")
        spec = importlib.util.spec_from_file_location("pop_params", path)
        _pop_params = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_pop_params)
    return _pop_params

def get_Bcur(Banc):
    params = get_pop_params()
    R = float(params.Nanc) / float(params.Ncur)
    Bcur = (Banc * (1.0 + (R - 1.0) * np.exp((-1.0 * params.time_of_change) / Banc))) / \
           (1 + (R - 1.0) * np.exp(-1.0 * params.time_of_change))
    return Bcur
