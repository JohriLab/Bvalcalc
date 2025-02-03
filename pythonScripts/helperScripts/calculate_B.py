import numpy as np
from constants import g, tract_len, r, u, Nanc, gamma_cutoff, h, t0, t1, t1half, t2, t3, t4, f0, f1, f2, f3


def calculate_exponent(t_start, t_end, U, a, b):
    """"
    Helper to calculate the exponent using "a" and "b"
    """
    E1 = ((U * a) 
            / ((1 - a) * (a - b) * (t_end - t_start))) * np.log((a + (t_end * (1 - a))) 
            / (a + (t_start * (1 - a))))
    E2 = -1.0 * ((U * b) 
            / ((1 - b) * (a - b) * (t_end - t_start))) * np.log((b + ((1 - b) * t_end)) 
            / (b + ((1 - b) * t_start)))
    return E1 + E2 # = E

def calculate_B(distance_to_element, length_of_element):
    """
    Calculate the B value for a single functional element at the focal site,
    summing over the DFE while consolidating the intermediate calculations.
    """
    # Calculate "a" and "b"
    C = (1.0 - np.exp(-2.0 * r * distance_to_element)) / 2.0
    U = length_of_element * 2 * u
    if g == 0:
        a = C
        b = C + r * length_of_element
    elif g > 0:
        threshold = distance_to_element + length_of_element < 0.5 * tract_len # Arbitrary threshold
        a = np.where(
            threshold, 
            C + (g * distance_to_element), #If TRUE
            C + g * tract_len #If FALSE
        )
        b = np.where(
            threshold,
            C + r * length_of_element + (g * (distance_to_element + length_of_element)), #If TRUE
            C + g * tract_len + r * length_of_element #If FALSE
        )

    # Calculate exponents for different time intervals
    E_f1 = calculate_exponent(t1half, t2, U, a, b)
    E_f2 = calculate_exponent(t2, t3, U, a, b)
    E_f3 = calculate_exponent(t3, t4, U, a, b)

    # Calculate E_bar (sum over DFE)
    E_bar = (
        f0 * 0.0
        + f1 * ((t1half - t1) / (t2 - t1)) * 0.0
        + f1 * ((t2 - t1half) / (t2 - t1)) * E_f1
        + f2 * E_f2
        + f3 * E_f3
    )

    # Return B
    return np.exp(-1.0 * E_bar)