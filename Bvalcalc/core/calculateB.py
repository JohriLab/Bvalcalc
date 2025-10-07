import numpy as np
from Bvalcalc.utils.dfe_helper import get_DFE_params
from scipy.optimize import root_scalar
from scipy.integrate import trapezoid

_params_cache: dict | None = None
_cache_args: tuple[str | None, bool, bool] | None = None

def get_params(
    params_path: str | None = None,
    gamma_dfe: bool = False,
    constant_dfe: bool = False,
):
    """
    Loads DFE parameters from the provided population genetic parameters file.
    Caches on (params_path, gamma_dfe, constant_dfe) and rebuilds whenever
    any of those three inputs change.
    """
    global _params_cache, _cache_args
    from Bvalcalc.utils.dfe_helper import GAMMA_DFE, CONSTANT_DFE
    
    # Get the actual params path that will be used
    actual_params_path = params_path
    if actual_params_path is None:
        import os
        actual_params_path = os.environ.get("BCALC_params")
    
    # Include global DFE state in cache key
    key = (actual_params_path, gamma_dfe or GAMMA_DFE, constant_dfe or CONSTANT_DFE)
    if _cache_args != key:
        _params_cache = get_DFE_params(params_path, gamma_dfe, constant_dfe)
        _cache_args = key
    return _params_cache

def calculateB_linear(distance_to_element: int, length_of_element: int, params: dict | None = None):
    """
    Calculate B due to purifying selection acting on a linked selected element of arbitrary length, assuming a constant crossover and gene conversion rate (analytical solution).
 
    Parameters 
    ----------
    distance_to_element: int
        Distance (bp) from the neutral site to the nearest edge of the selected element.
    length_of_element: int
        Length (bp) of the selected element.
    params : dict
        Required parameters from ``get_params()``, only kept as default (None) when being called by CLI,
        in which case parameters are sourced from the params file directly.
    """    
    with np.errstate(divide='ignore', invalid='ignore'):
        if params is None:
            params = get_params()
        r, u, g, k, t1, t1half, t2, t3, t4, f1, f2, f3, f0, t_constant = params["r"], params["u"], params["g"], params["k"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f1"], params["f2"], params["f3"], params["f0"], params["t_constant"]

        C = (1.0 - np.exp(-2.0 * r * distance_to_element)) / 2.0 # cM
        U = length_of_element * u
        if g == 0:
            a = C # RECOMBINATION IN Y
            b = C + (r * length_of_element) # RECOMBINATION IN X
        elif g > 0:
            a, b = get_a_b_with_GC(C, distance_to_element, length_of_element, params)

        if t_constant: #If --constant_dfe is active
            E_constant = calculate_exponent(t_constant, t_constant, U, a, b)
            B = np.exp(-1.0 * E_constant)
            return np.where(length_of_element == 0, 1.0, B)
        
        E_f1 = calculate_exponent(t1half, t2, U, a, b)
        E_f2 = calculate_exponent(t2, t3, U, a, b)
        E_f3 = calculate_exponent(t3, t4, U, a, b)

        E_bar = ( # Sum over the DFE
            f0 * 0.0
            + f1 * ((t1half - t1) / (t2 - t1)) * 0.0
            + f1 * ((t2 - t1half) / (t2 - t1)) * E_f1
            + f2 * E_f2
            + f3 * E_f3)

        B = np.exp(-1.0 * E_bar)
        
    return np.where(length_of_element == 0, 1.0, B)

def calculateB_recmap(distance_to_element, length_of_element, 
                      rec_distances = None, rec_lengths = None, 
                      gc_distances = None, gc_lengths = None, params = None):
    """
    Calculate the B value WITH REC MAP for a single functional element at the focal site,
    summing over the DFE while consolidating the intermediate calculations.
    """    
    with np.errstate(divide='ignore', invalid='ignore'):
        if params is None:
            params = get_params()
        r, u, g, k, t1, t1half, t2, t3, t4, f1, f2, f3, f0, t_constant = params["r"], params["u"], params["g"], params["k"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f1"], params["f2"], params["f3"], params["f0"], params["t_constant"]
        # rec_distances is the length of the element * rec rate in each spanned region. 
        
        if rec_distances is not None:
            rec_adjusted_length_of_element = rec_lengths 
            rec_adjusted_distance_to_element = rec_distances
        else:
            rec_adjusted_length_of_element = length_of_element
            rec_adjusted_distance_to_element = distance_to_element
        
        if gc_distances is not None:
            local_g = (gc_lengths + gc_distances)/(length_of_element + distance_to_element) * g
        else:
            local_g = g
        
        C = (1.0 - np.exp(-2.0 * r * rec_adjusted_distance_to_element)) / 2.0 # cM
        U = length_of_element * u
        if g == 0:
            a = C
            b = C + r * rec_adjusted_length_of_element # cM
        elif g > 0:
             a, b = get_a_b_with_GC_andMaps(C, y=distance_to_element, l=length_of_element, 
                                            rec_l=rec_adjusted_length_of_element, local_g = local_g, params=params)

        if t_constant: #If --constant_dfe is active
            E_constant = calculate_exponent(t_constant, t_constant, U, a, b)
            B = np.exp(-1.0 * E_constant)
            return np.where(length_of_element == 0, 1.0, B)
        
        E_f1 = calculate_exponent(t1half, t2, U, a, b)
        E_f2 = calculate_exponent(t2, t3, U, a, b)
        E_f3 = calculate_exponent(t3, t4, U, a, b)

        E_bar = ( # Sum over the DFE
            f0 * 0.0
            + f1 * ((t1half - t1) / (t2 - t1)) * 0.0
            + f1 * ((t2 - t1half) / (t2 - t1)) * E_f1
            + f2 * E_f2
            + f3 * E_f3)    

        B = np.exp(-1.0 * E_bar)
        
    return np.where(length_of_element == 0, 1.0, B)

def calculateB_unlinked(unlinked_L: int, params: dict | None = None):
    """
    Calculate B due to purifying selection at unlinked sites.

    Parameters
    ----------
    unlinked_L : float
        Cumulative count of selected sites in unlinked regions.
    params : dict
        Required parameters from ``get_params()``, only kept as default (None) when being called by CLI,
        in which case parameters are sourced from the params file directly.
    """
    if params is None:
        params = get_params()

    u, t1, t1half, t2, t3, t4, f0, f1, f2, f3, t_constant = params["u"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f0"], params["f1"], params["f2"], params["f3"], params["t_constant"]
    
    if t_constant: #If --constant_dfe is active    

        unlinked_B  = np.exp(-8 * u * 1.0 * unlinked_L * (t_constant/(1 + t_constant)**2))
        # unlinked_B  = np.exp(-8 * u * 1.0 * unlinked_L * t_constant) ## THIS IS EQ. XX APPROXIMATION IN THE MANUSCRIPT
        return unlinked_B    

    f1_above_cutoff = f1 * ((t1half - t1) / (t2 - t1))

    sum_f1 = (f1_above_cutoff / (t2 - t1half)) * (np.log((1 + t2) /(1 + t1half)) + (1 / (1 + t2)) - (1 / (1 + t1half)))
    sum_f2 = (f2 / (t3 - t2)) * (np.log((1 + t3) /(1 + t2)) + (1 / (1 + t3)) - (1 / (1 + t2)))
    sum_f3 = (f3 / (t4 - t3)) * (np.log((1 + t4) /(1 + t3)) + (1 / (1 + t4)) - (1 / (1 + t3)))
    
    unlinked_B  = np.exp(-8 * u * 1.0 * unlinked_L * (sum_f1 + sum_f2 + sum_f3))

    return unlinked_B


##



## Helper functions

def calculate_exponent(t_start, t_end, U, a, b):
    """"
    Helper to calculate the exponent using "a" and "b"
    """
    a, b, U = np.asarray(a), np.asarray(b), np.asarray(U)

    if U.size == 0: return 0 # If e.g. f1 proportion is 0, no need to calculate exponent
    
    if t_end == t_start: # If --constant_dfe
        E = (U / (a - b)) * (
            a / (a + (1 - a) * t_start) -
            b / (b + (1 - b) * t_start)
        )
    else: # Using discretized DFE (f0,f1,f2,f3 or --gamma_dfe)
        E1 = ((U * a) 
                / ((1 - a) * (a - b) * (t_end - t_start))) * np.log((a + (t_end * (1 - a))) 
                / (a + (t_start * (1 - a))))
        E2 = -1.0 * ((U * b) 
                / ((1 - b) * (a - b) * (t_end - t_start))) * np.log((b + ((1 - b) * t_end)) 
                / (b + ((1 - b) * t_start)))
    
        E = np.asarray(E1 + E2)

    rec_0_mask = np.isclose(a, b)  # Get mask for where recombination rate = 0 within the gene
    if rec_0_mask.any(): # 4a) If a_arr is scalar (0‐d), compute limit once as scalar
        if a.ndim == 0:
            limit_factor = (1 / ((t_end - t_start)*(1-a)**2)) * ( # Calculate exponent with 0 recombination between gene and site, avoiding limits
                np.log((a + (1 - a) * t_end) 
                       / (a + (1 - a) * t_start))
                + a / (a + (1 - a) * t_end)
                - a / (a + (1 - a) * t_start))
            if t_start == t_end: 
                limit_factor = t_start / (a + (1 - a) * t_start)**2 # If --constant_dfe
                return U * limit_factor # E is numpy scalar when t_constant and a = b
            # Broadcast scalar limit_factor to all masked positions
            E[rec_0_mask] = U[rec_0_mask] * limit_factor  # Get corresponding U for the numerator and plug back into E array to replace nan's

        else: # 4b) If a_arr is array, compute limit for each masked element
            ae = a[rec_0_mask]  # array of a_i where a_i ≈ b_i
            limit_factor = (1 / ((t_end - t_start)*(1-ae)**2)) * ( # Calculate exponent with 0 recombination between gene and site, avoiding limits
                np.log((ae + (1 - ae) * t_end) 
                       / (ae + (1 - ae) * t_start))
                + ae / (ae + (1 - ae) * t_end)
                - ae / (ae + (1 - ae) * t_start))
            if t_start == t_end: limit_factor = t_start / (ae + (1 - ae) * t_start)**2 # If --constant_dfe
            ## REPLACED BELOW WITH THE NEW LINE TO FIX FAR GENE ISSUE, MAY NEED TO REVERT
            E[rec_0_mask] = U[rec_0_mask] * limit_factor
            # Match array of limit_factor to corresponding positions in E (where rec_0_mask has True);'l;'l''
            # if len(rec_0_mask[False]) == 0:
            #     # print(f"Need to fix --gene when r = 0, see calculateB ~line 176") Fixed??
            #     E[rec_0_mask] = U * limit_factor
            # else:
            #     E[rec_0_mask] = U[rec_0_mask] * limit_factor  # Get corresponding U for the numerator and plug back into E array to replace nan's

    return E

def get_a_b_with_GC(C, y, l, params=None):
        with np.errstate(divide='ignore', invalid='ignore'):
            if params is None:
                params = get_params()
            r, u, g, k, t1, t1half, t2, t3, t4, f1, f2, f3, f0 = params["r"], params["u"], params["g"], params["k"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f1"], params["f2"], params["f3"], params["f0"]
            proportion_nogc_a = np.where(k < y + l, # When GC includes neutral site, this is proportion of the gene it includes
                                        np.maximum((0.5*(k-y)/l), 0),
                                        1-(y + l)/(2 * k)
                                        )

            proportion_nogc_b = np.where(k < y + l, # When GC includes gene site, this is probability the tract includes neutral site of interest 
                                    1/(2*k) * np.maximum(k-y+1,0) * np.maximum(k - y, 0) / l, # When overshooting not possible
                                    (k - y - 0.5 * l) / k) # When overshooting possible
            
        
        a = np.where(k < y, 
            C + (2 * g * k), # Probability of GC on neutral site, where overlap with element not possible
            C + (2 * g * (y) + # When overlap possible this is probability gc is in neutral but doesn't include any of element
                g * (k - y) * # Probability gc is in neutral and includes some element (remaining probability from above)
                (1 - proportion_nogc_a) # Proportion of gene that gc breaks linkage with when it includes some element
        ))
        b = C + (r * l) + (2 * g * k) * (1 - (1-proportion_nogc_a)*proportion_nogc_b) #* prop k out

        return a, b

def get_a_b_with_GC_andMaps(C, y, l, rec_l, local_g, params=None):
        if params is None:
            params = get_params()
        r, u, g, k, t1, t1half, t2, t3, t4, f1, f2, f3, f0 = params["r"], params["u"], params["g"], params["k"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f1"], params["f2"], params["f3"], params["f0"]
        with np.errstate(divide='ignore', invalid='ignore'):
            proportion_nogc_a = np.where(k < y + l, # When GC includes neutral site, this is proportion of the gene it includes
                                        np.maximum((0.5*(k-y)/l), 0),
                                        ((y) * (2 * k - (y + l)))/(2 * k * y)
                                        )

            proportion_nogc_b = np.where(k < y + l, # When GC includes gene site, this is probability the tract includes neutral site of interest 
                                    1/(2*k) * np.maximum(k-y+1,0) * np.maximum(k - y, 0) / l,
                                    (k - y - 0.5 * l) / k)
        
        a = np.where(k < y, 
            C + (2 * local_g * k), # Probability of GC on neutral site, where overlap with element not possible
            C + (2 * local_g * (y) + # When overlap possible this is probability gc is in neutral but doesn't include any of element
                local_g * (k - y) * # Probability gc is in neutral and includes some element (remaining probability from above)
                (1 - proportion_nogc_a) # Proportion of gene that gc breaks linkage with when it includes some element
        ))
        b = C + (r * rec_l) + (2 * local_g * k) * (1 - (1-proportion_nogc_a)*proportion_nogc_b) #* prop k out

        return a, b

def calculateB_hri(distant_B, interfering_L, params: dict | None = None):
    """
    Calculate B' (B accounting for Hill-Robertson interference effects) for a non-recombining region containing selected sites. 

    This is a diploid implementation of Eq. 12 from Becher and Charlesworth (2025), see the relevant supplement in the B-value calculator manuscript.

    Parameters
    ----------
    distant_B : float or array-like
        The background B value in the HRI region, from less-linked selected sites outside the interference region (i.e. elsewhere on the chromosome, other chromosomes).
    interfering_L : float or array-like
        The cumulative length (bp) of interfering selected sites in the HRI region.
    params : dict
        Required parameters from ``get_params()``, only kept as default (None) when being called by CLI,
        in which case parameters are sourced from the params file directly.
    """
    if params is None:
        params = get_params()

    Nanc, u, f1, f2, t_constant = params["Nanc"], params["u"], params["f1"], params["f2"], params["t_constant"]

    distant_B = np.atleast_1d(distant_B).astype(float)
    interfering_L = np.atleast_1d(interfering_L).astype(float)

    scalar_input = distant_B.shape == () or distant_B.shape == (1,)

    # Early return: if no interfering L in inputs, B' = distant_B
    if np.all(interfering_L == 0):
        return distant_B[0] if scalar_input else distant_B

    N0 = distant_B * Nanc
    h = 0.5
    u = 2 * u
    u1 = f1 * u
    u2 = f2 * u
    u_total = u1 + u2

    E_X2_f1 = (1**2 + 1*10 + 10**2) / 3
    E_X2_f2 = (10**2 + 10*100 + 100**2) / 3

    t_sq1 = (h**2 * E_X2_f1) / (4 * N0**2)
    t_sq2 = (h**2 * E_X2_f2) / (4 * N0**2)
    t = np.sqrt((u1 * t_sq1 + u2 * t_sq2) / u_total)

    if t_constant:
        t = t_constant
        u_total = u

    gamma = 2 * N0 * t
    U = u_total * interfering_L
    alpha2 = 2 * N0 * U
    kappa = 1.0

    def eq4(B, U, gamma, t):
        exp_term = np.exp(-gamma * B)
        num = 0.5 * U * (1 - exp_term)**3
        denom = t * (1 + kappa * exp_term)**3
        return -np.log(B) - num / denom

    def solve_eq4_batched(U, gamma, t, n=500):
        Bgrid = np.linspace(1e-10, 1.0, n)[None, :]
        U = np.asarray(U).reshape(-1, 1)
        gamma = np.asarray(gamma).reshape(-1, 1)
        t = np.asarray(t).reshape(-1, 1)

        fvals = eq4(Bgrid, U, gamma, t)
        signs = np.sign(fvals)
        crossing = np.diff(signs, axis=1) < 0
        idx = np.argmax(crossing, axis=1)

        B_left = Bgrid[0, idx]
        B_right = Bgrid[0, idx + 1]
        f_left = fvals[np.arange(len(U)), idx]
        f_right = fvals[np.arange(len(U)), idx + 1]

        B_root = B_left - f_left * (B_right - B_left) / (f_right - f_left)
        return B_root

    Bval = solve_eq4_batched(U, gamma, t)

    def eq5_vectorized(B, alpha2, gamma, Tmax=100.0, n_steps=2000):
        x = np.linspace(0, Tmax, n_steps)[None, :]  # shape (1, n_steps)
        dx = x[0, 1] - x[0, 0]

        B = B[:, None]
        alpha2 = alpha2[:, None]
        gamma = gamma[:, None]

        f1 = 1 - np.exp(-gamma * B)
        f2 = 1 + kappa * np.exp(-gamma * B)
        A = f1 / f2
        c = 0.5 * alpha2 / gamma * A**3
        d = 2 * gamma * B * (f2 / f1)

        x_broadcasted = np.broadcast_to(x, (B.shape[0], x.shape[1]))
        gx = np.exp(c * (1 - np.exp(-d * x_broadcasted))**2)
        cumI = np.cumsum((gx[:, :-1] + gx[:, 1:]) * 0.5 * dx, axis=1)
        cumI = np.hstack([np.zeros((gx.shape[0], 1)), cumI])

        hx = np.exp(-B * cumI)
        Bprime = B[:, 0] * trapezoid(hx, x[0], axis=1)
        return Bprime

    Bprime = eq5_vectorized(Bval, alpha2, gamma)
    return Bprime[0] if scalar_input else Bprime