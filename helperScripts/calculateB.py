import numpy as np
from ExampleParams import g, k, r, u, t1, t1half, t2, t3, t4, f0, f1, f2, f3

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

def get_a_b_with_GC(C, y, l):
        proportion_nogc_b = np.where(k < y + l, # When GC includes gene site, this is probability the tract includes neutral site of interest 
                                   1/(2*k) * np.maximum(k-y+1,0) * np.maximum(k - y, 0) / l,
                                   (k - y - 0.5 * l) / k)
        
        proportion_nogc_a = np.where(k < y + l, # When GC includes neutral site, this is proportion of the gene it includes
                                    np.maximum((0.5*(k-y)/l), 0),
                                    ((y) * (2 * k - (y + l)))/(2 * k * y)
                                    )

        a = np.where(k < y, 
            C + (2 * g * k), # Probabiliity of GC on neutral site, where overlap with element not possible
            C + (2 * g * (y) + # When overlap possible this is probability gc is in neutral but doesn't include any of element
                g * (k - y) * # Probability gc is in neutral and includes some element (remaining probability from above)
                (1 - proportion_nogc_a) # Proportion of gene that gc breaks linkage with when it includes some element
        ))
        b = C + (r * l) + (2 * g * k) * (1 -  proportion_nogc_b) #* prop k out

        return a, b

def calculateB_linear(distance_to_element, length_of_element):
    """
    Calculate the B value for a single functional element at the focal site,
    summing over the DFE while consolidating the intermediate calculations.
    """    
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

    return np.exp(-1.0 * E_bar) # Return B


def calculateB_recmap(distance_to_element, length_of_element, 
                      rec_distances = None, rec_lengths = None, 
                      gc_distances = None, gc_lengths = None):
    """
    Calculate the B value WITH REC MAP for a single functional element at the focal site,
    summing over the DFE while consolidating the intermediate calculations.
    """    
    # rec_distances is the length of the element * rec rate in each spanned region. 
    if rec_distances is not None:
        rec_adjusted_length_of_element = rec_lengths 
        rec_adjusted_distance_to_element = rec_distances
    else:
        rec_adjusted_length_of_element = length_of_element
        rec_adjusted_distance_to_element = distance_to_element
    
    if gc_distances is not None:
        gc_adjusted_length_of_element = gc_lengths 
        gc_adjusted_distance_to_element = gc_distances
        g_k = (gc_lengths / length_of_element) * g * k
    else:
        gc_adjusted_length_of_element = length_of_element
        gc_adjusted_distance_to_element = distance_to_element
        g_k = g * k
        
    C = (1.0 - np.exp(-2.0 * r * rec_adjusted_distance_to_element)) / 2.0 # cM
    U = length_of_element * u
    if g == 0:
        a = C
        b = C + r * rec_adjusted_length_of_element # cM
    elif g > 0:
        threshold = gc_adjusted_distance_to_element + gc_adjusted_length_of_element < 0.5 * k # Arbitrary threshold
        a = np.where(
            threshold, 
            C + (g * gc_adjusted_distance_to_element), #If TRUE
            C + g_k #If FALSE
        )
        b = np.where(
            threshold,
            C + r * rec_adjusted_length_of_element + (g * (gc_adjusted_distance_to_element + gc_adjusted_length_of_element)), #If TRUE
            C + g_k + r * rec_adjusted_length_of_element #If FALSE
        )

    E_f1 = calculate_exponent(t1half, t2, U, a, b)
    E_f2 = calculate_exponent(t2, t3, U, a, b)
    E_f3 = calculate_exponent(t3, t4, U, a, b)

    E_bar = ( # Sum over the DFE
        f0 * 0.0
        + f1 * ((t1half - t1) / (t2 - t1)) * 0.0
        + f1 * ((t2 - t1half) / (t2 - t1)) * E_f1
        + f2 * E_f2
        + f3 * E_f3)

    return np.exp(-1.0 * E_bar) # Return B

def calculateB_linear_oldGC(distance_to_element, length_of_element):
    """
        Not currently used in main Bvalcalc function, 
        Here for reference, Parul's gene conversion equation.
    """    
    C = (1.0 - np.exp(-2.0 * r * distance_to_element)) / 2.0 # cM
    U = length_of_element * u
    if g == 0:
        a = C
        b = C + (r * length_of_element) # cM
    elif g > 0:
        threshold = distance_to_element + length_of_element < 0.5 * k # Arbitrary threshold
        a = np.where(
            threshold, 
            C + (g * distance_to_element), #pGC happens outside the element to BREAK linkage, lowers because sometimes both GC
            C + g * k #pGC happens outside element to BREAK linkage
            ) 
        b = np.where(
            threshold,
            C + (r * length_of_element) + (g * (distance_to_element + length_of_element)), #pGC happens within element to BREAK linkage, IS IT LOWER???
            C + (r * length_of_element) + (g * k) #pGC happens within element to BREAK linkage,
        )

    E_f1 = calculate_exponent(t1half, t2, U, a, b)
    E_f2 = calculate_exponent(t2, t3, U, a, b)
    E_f3 = calculate_exponent(t3, t4, U, a, b)

    E_bar = ( # Sum over the DFE
        f0 * 0.0
        + f1 * ((t1half - t1) / (t2 - t1)) * 0.0
        + f1 * ((t2 - t1half) / (t2 - t1)) * E_f1
        + f2 * E_f2
        + f3 * E_f3)

    return np.exp(-1.0 * E_bar) # Return B