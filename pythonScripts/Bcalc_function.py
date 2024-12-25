import math
import numpy as np
from constants import g, tract_len, r, u, l, U, Nanc, gamma_cutoff, h, t0, t1, t1half, t2, t3, t4, f0, f1, f2, f3

# Function to calculate "a" and "b"
def calculate_a_and_b(distance_to_element, length_of_element):
    C = (1.0 - math.exp(-2.0 * r * distance_to_element)) / 2.0
    if g == 0:
        a = C
        b = C + r * length_of_element
    elif g > 0:
        if distance_to_element + length_of_element < 0.5 * tract_len:  # Arbitrary threshold
            a = C + (g * distance_to_element)
            b = C + r * length_of_element + (g * (distance_to_element + length_of_element))
        else:
            a = g * tract_len + C
            b = g * tract_len + r * length_of_element + C
    return a, b

# Function to calculate the exponent using "a" and "b"
def calculate_exponent(t_start, t_end, distance_to_element, length_of_element):
    a, b = calculate_a_and_b(distance_to_element, length_of_element)
    E1 = ((U * a) / ((1 - a) * (a - b) * (t_end - t_start))) * math.log((a + (t_end * (1 - a))) / (a + (t_start * (1 - a))))
    E2 = -1.0 * ((U * b) / ((1 - b) * (a - b) * (t_end - t_start))) * math.log((b + ((1 - b) * t_end)) / (b + ((1 - b) * t_start)))
    return E1 + E2

# Function to calculate B for a single functional element
def calculate_B(distance_to_element, length_of_element):
    E_f1 = calculate_exponent(t1half, t2, distance_to_element, length_of_element)
    E_f2 = calculate_exponent(t2, t3, distance_to_element, length_of_element)
    E_f3 = calculate_exponent(t3, t4, distance_to_element, length_of_element)
    E_bar = (
        f0 * 0.0
        + f1 * ((t1half - t1half) / (t2 - t1half)) * 0.0
        + f1 * ((t2 - t1half) / (t2 - t1half)) * E_f1
        + f2 * E_f2
        + f3 * E_f3
    )
    return math.exp(-1.0 * E_bar)

# Vectorized version of calculate_B
vectorized_B = np.vectorize(calculate_B)
