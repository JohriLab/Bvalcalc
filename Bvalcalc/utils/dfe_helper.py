import scipy.stats as st
import os
import importlib.util
import numpy as np
from typing import Dict, Any

GAMMA_DFE = False  # Default, instead of prop injected
CONSTANT_DFE = False # Default, instead of prop injected
CUSTOM_DFE = False # Default, instead of prop injected

def get_DFE_params(params_path: str | None = None, gamma_dfe: bool = False, constant_dfe: bool = False, custom_dfe: bool = False) -> Dict[str, Any]:
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

    # 5. Set derived parameters and thresholds
    params["gamma_cutoff"] = 5
    s_cutoff = params["gamma_cutoff"] / (2.0 * params["Nanc"]) # Convert 2Ns cutoff to s cutoff
    params["t0"] = 0.0
    Nanc = params["Nanc"]
    h = params["h"]
    # Calculate generation-scale thresholds
    params["t1"] = h * (1.0 / (2.0 * Nanc))
    params["t1half"] = h * (params["gamma_cutoff"] / (2.0 * Nanc))
    params["t2"] = h * (10.0 / (2.0 * Nanc))
    params["t3"] = h * (100.0 / (2.0 * Nanc))
    params["t4"] = h * 1.0



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

        f0, f_x, s_edges = gammaDFE_to_discretized(
            mean,
            shape,
            prop_syn,
            s_cutoff=s_cutoff,
        )

        params["f0"] = f0
        params["f_x"] = f_x
        params["t_edges"] = h * s_edges



    elif CUSTOM_DFE or custom_dfe is not False: # The CUSTOM_DFE is prop injected by CLI, custom_dfe is provided by API
        s_breaks = getattr(pop, 's_breaks', None)
        bin_props = getattr(pop, 'bin_proportions', None)
        if len(s_breaks) != len(bin_props) + 1:
            raise ValueError(
                "Length of s_breaks must be one more than length of bin_proportions when --custom_dfe is active. Breaks define the edges of the bins, so there must be one more break than bin."
            )
        if s_breaks[0] != 0 or s_breaks[-1] != 1:
            raise ValueError(
                "s_breaks must start with 0 and end with 1 when --custom_dfe is active. This defines the full range of selection coefficients from neutral (0) to fully homozygous lethal (1). If your distribution does not include neutral or fully lethal mutations, set the edge proportions to 0"
            )
        if any((x < 0 or x > 1) for x in s_breaks):
            raise ValueError(
                "s_breaks defines deleterious DFE values, for simplicity here, positive values provided to s_breaks represent the strength of purifying selection (a value of 1 is homozygous lethal, equal to s = -1)."
            )
        if not np.isclose(sum(bin_props), 1.0):
            raise ValueError(
                "bin_proportions must sum to 1 when --custom_dfe is active. These define the proportion of mutations in each bin defined by s_breaks."
            )
        
        ## Here we collapse all proportions below s_cutoff to a single f0 proportion.
        from .dfe_helper import customDFE_getf0
        f0, s_breaks, bin_props = customDFE_getf0(s_breaks, bin_props, s_cutoff)

        params["f0"] = f0
        params["f_x"] = bin_props
        params["t_edges"] = h * s_breaks
    
    else:
        from .dfe_helper import legacyDFE_to_bins
        params = legacyDFE_to_bins(params)

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


def gammaDFE_to_discretized(mean, shape, proportion_synonymous,
                            s_cutoff=None):

    if mean <= 0 or shape <= 0:
        raise ValueError("`mean` and `shape` must be positive.")

    if not (0 <= proportion_synonymous < 1):
        raise ValueError("`proportion_synonymous` must be in [0, 1).")

    theta = mean / shape
    dist = st.gamma(a=shape, scale=theta)

    # ------------------------------------------------------------------
    # fixed log-spaced bins in `s` space from 0 (neutral) to 1 (homozygous lethal)
    # ------------------------------------------------------------------

    s_edges = np.array([
        0.0,
        1e-8,
        1e-7,
        1e-6,
        1e-5,
        1e-4,
        1e-3,
        1e-2,
        1e-1,
        1.0,
    ], dtype=float)

    # discretized proportions
    f_x = np.zeros(len(s_edges) - 1, dtype=float)

    for i in range(len(f_x)):

        left = s_edges[i]
        right = s_edges[i + 1]

        f_x[i] = dist.cdf(right) - dist.cdf(left)

    # fold tail mass (>1) into final bin
    f_x[-1] += 1.0 - dist.cdf(s_edges[-1])

    # scale to nonsynonymous fraction
    scale = 1.0 - proportion_synonymous
    f_x *= scale

    # collapse effectively-neutral mass into f0
    f0 = proportion_synonymous

    if s_cutoff is not None:
        gamma_f0 = dist.cdf(s_cutoff) * scale
        f0 += gamma_f0

        for i, prop in enumerate(f_x):
            left = s_edges[i]
            right = s_edges[i + 1]

            # Entire bin below cutoff
            if right <= s_cutoff:
                f_x[i] = 0.0

            # Bin crosses cutoff
            elif left < s_cutoff < right:
                above_fraction = ((right - s_cutoff) / (right - left))
                f_x[i] *= above_fraction

        # renormalize remaining selected mass
        remaining = f_x.sum()

        if remaining > 0:
            f_x *= (scale - gamma_f0) / remaining

    print(f"Converting gamma distribution to discretized DFE")
    print(f"Gamma params: mean s = {mean:.8f}, shape = {shape}, scale = {theta:.8f}")
    print(f"s_edges, selection coefficient breakpoint for each bin = {s_edges}")
    print(f"s_cutoff, minimum selection coefficient for which BGS is calculated (below added to f0) = {s_cutoff:.6f}")
    print(f"f0, effectively neutral proportion in selected region = {f0:.6f}, of which {gamma_f0 if s_cutoff is not None else 0.0:.6f} is from the gamma distribution below cutoff")
    print(f"f_x, proportion in each remaining bin, excluding f0 = {np.array2string(f_x, formatter={'float_kind': lambda x: f'{x:.6f}'})}")

    return (f0, np.array(f_x, dtype=float), np.array(s_edges, dtype=float))

def customDFE_getf0(s_breaks, bin_props, s_cutoff):
    """
    Collapse all DFE mass below s_cutoff into f0, while preserving
    the remaining bins above the cutoff.

    Returns
    -------
    f0 : float
        Total proportion below cutoff.
    new_breaks : np.ndarray
        Bin edges from s_cutoff upward.
    new_props : np.ndarray
        Bin proportions corresponding to new_breaks.
    """

    s_breaks = np.array(s_breaks, dtype=float)
    bin_props = np.array(bin_props, dtype=float)

    new_breaks = [s_cutoff]
    new_props = []

    f0 = 0.0

    for i, prop in enumerate(bin_props):

        left = s_breaks[i]
        right = s_breaks[i + 1]

        # Entire bin below cutoff
        if right <= s_cutoff:
            f0 += prop

        # Bin crosses cutoff
        elif left < s_cutoff < right:

            below_fraction = (s_cutoff - left) / (right - left)
            above_fraction = 1.0 - below_fraction

            f0 += prop * below_fraction
            new_props.append(prop * above_fraction)

            new_breaks.append(right)

        # Entire bin above cutoff
        else:
            new_props.append(prop)
            new_breaks.append(right)
    
    print(f"Custom DFE parameters active (overwriting f1-f3)")
    print(f"s_edges, homozygous selection coefficient breakpoint for each bin =", np.array2string(np.array(new_breaks, dtype=float),precision=3))
    print(f"s_cutoff, minimum selection coefficient for which BGS is calculated (below added to f0) = {s_cutoff:.6f}")
    print(f"f0, effectively neutral proportion in selected region = {f0:.6f} which is not used in BGS calculations")
    print(f"f_x, proportion in each remaining bin, excluding f0 =", np.array2string(np.array(new_props, dtype=float),formatter={'float_kind': lambda x: f'{x:.6f}'}))

    return ( # Return s clean arrays for exporting to calcB
        f0,
        np.array(new_breaks, dtype=float),
        np.array(new_props, dtype=float),
    )


def legacyDFE_to_bins(params):
    """
    Convert legacy f1/f2/f3 discretized DFE into
    unified f_x / t_edges representation.
    """

    t1 = params["t1"]
    t1half = params["t1half"]
    t2 = params["t2"]
    t3 = params["t3"]
    t4 = params["t4"]

    f1 = params["f1"]
    f2 = params["f2"]
    f3 = params["f3"]

    # remove effectively-neutral portion of f1
    f1_selected = f1 * ((t2 - t1half) / (t2 - t1))

    params["f_x"] = np.array([
        f1_selected,
        f2,
        f3
    ], dtype=float)

    params["t_edges"] = np.array([
        t1half,
        t2,
        t3,
        t4
    ], dtype=float)

    print(f"Basic DFE parameters active by default (see --constant_dfe, --gamma_dfe, --custom_dfe)")

    return params