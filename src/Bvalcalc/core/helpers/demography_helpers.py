import os
import numpy as np
import importlib.util
from typing import TYPE_CHECKING

# Load and cache the parameters
if TYPE_CHECKING:
    Nanc: float; Ncur: float; time_of_change: float
spec = importlib.util.spec_from_file_location("pop_params", os.environ["BCALC_POP_PARAMS"])
_pop = importlib.util.module_from_spec(spec); spec.loader.exec_module(_pop)
for v in ['Nanc', 'Ncur', 'time_of_change']: globals()[v] = getattr(_pop, v)

def get_Bcur(Banc):
    R = Nanc / Ncur
    Bcur = (Banc * (1 + (R - 1) * np.exp(-(time_of_change/(2*Ncur)) / Banc))) / \
           (1 + (R - 1) * np.exp(-(time_of_change/(2*Ncur)))) # Denominator
    return Bcur
