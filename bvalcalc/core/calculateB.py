import numpy as np
from bvalcalc.utils.dfe_helper import get_DFE_params
from scipy.integrate import quad
import sys

_params_cache: dict | None = None

def get_params(params_path: str | None = None):
    """
    Loads DFE parameters from the provided population genetic parameters file.

    Parameters 
    ----------
    params_path: str
        Path to Params.py file
    """
    global _params_cache
    if _params_cache is None:
        _params_cache = get_DFE_params(params_path)
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
        r, u, g, k, t1, t1half, t2, t3, t4, f1, f2, f3, f0 = params["r"], params["u"], params["g"], params["k"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f1"], params["f2"], params["f3"], params["f0"]
        C = (1.0 - np.exp(-2.0 * r * distance_to_element)) / 2.0 # cM
        U = length_of_element * u
        if g == 0:
            a = C # RECOMBINATION IN Y
            b = C + (r * length_of_element) # RECOMBINATION IN X
        elif g > 0:
            a, b = get_a_b_with_GC(C, distance_to_element, length_of_element)

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
        r, u, g, k, t1, t1half, t2, t3, t4, f1, f2, f3, f0 = params["r"], params["u"], params["g"], params["k"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f1"], params["f2"], params["f3"], params["f0"]
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
                                            rec_l=rec_adjusted_length_of_element, local_g = local_g)

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
    Calculate B due to purifying selection at unlinked sites (numerical integration over DFE).

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

    u, t1, t1half, t2, t3, t4, f0, f1, f2, f3 = params["u"], params["t1"], params["t1half"], params["t2"], params["t3"], params["t4"], params["f0"], params["f1"], params["f2"], params["f3"]

    def integrate_unlinked_B(f, t_a, t_b):
        """
        Gets B for uniform DFE from t_a to t_b (e.g. t_a=0.01 to t_b=0.5) with frequency f 
        For example, 30% of the DFE is f3 so f = 0.3, and t_a = 0.0025 and t_b = 0.5 (s = -1, h = 0.5). 
        """
        integrand = lambda t: np.exp(-8 * u * unlinked_L * f3 * t / (1 + t)**2)
        integral_f, _ = quad(integrand, t_a, t_b)
        B_f = integral_f / (t_b - t_a)

        return B_f

    B_from_f1 = integrate_unlinked_B(f1, t1half, t2)
    B_from_f2 = integrate_unlinked_B(f2, t2, t3)
    B_from_f3 = integrate_unlinked_B(f3, t3, t4)

    return B_from_f1 * B_from_f2 * B_from_f3


## Helper functions

def calculate_exponent(t_start, t_end, U, a, b):
    """"
    Helper to calculate the exponent using "a" and "b"
    """
    if a == b: # If recombination rate = 0 between the gene and the distance from the gene
        E = U / (t_end - t_start) * (np.log((a + (1 - a) * t_end) / 
                                    (a + (1 - a) * t_start)) / 
                                    (1 - a) ** 2 + 
        (a / (1 - a)) * ((1 - t_end) / 
                     (a + (1 - a) * t_end) - 
                    (1 - t_start) / 
                    (a + (1 - a) * t_start)))
        return E
    else:
        E1 = ((U * a) 
                / ((1 - a) * (a - b) * (t_end - t_start))) * np.log((a + (t_end * (1 - a))) 
                / (a + (t_start * (1 - a))))
        E2 = -1.0 * ((U * b) 
                / ((1 - b) * (a - b) * (t_end - t_start))) * np.log((b + ((1 - b) * t_end)) 
                / (b + ((1 - b) * t_start)))
        return E1 + E2 # = E

def get_a_b_with_GC(C, y, l):
        with np.errstate(divide='ignore', invalid='ignore'):
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

def get_a_b_with_GC_andMaps(C, y, l, rec_l, local_g):
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