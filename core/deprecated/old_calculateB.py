# def calculateB_linear_oldGC(distance_to_element, length_of_element):
#     """
#         Not currently used in main Bvalcalc function, 
#         Here for reference, Parul's gene conversion equation.
#     """    
#     C = (1.0 - np.exp(-2.0 * r * distance_to_element)) / 2.0 # cM
#     U = length_of_element * u
#     if g == 0:
#         a = C
#         b = C + (r * length_of_element) # cM
#     elif g > 0:
#         threshold = distance_to_element + length_of_element < 0.5 * k # Arbitrary threshold
#         a = np.where(
#             threshold, 
#             C + (g * distance_to_element), #pGC happens outside the element to BREAK linkage, lowers because sometimes both GC
#             C + g * k #pGC happens outside element to BREAK linkage
#             ) 
#         b = np.where(
#             threshold,
#             C + (r * length_of_element) + (g * (distance_to_element + length_of_element)), #pGC happens within element to BREAK linkage, IS IT LOWER???
#             C + (r * length_of_element) + (g * k) #pGC happens within element to BREAK linkage,
#         )

#     E_f1 = calculate_exponent(t1half, t2, U, a, b)
#     E_f2 = calculate_exponent(t2, t3, U, a, b)
#     E_f3 = calculate_exponent(t3, t4, U, a, b)

#     E_bar = ( # Sum over the DFE
#         f0 * 0.0
#         + f1 * ((t1half - t1) / (t2 - t1)) * 0.0
#         + f1 * ((t2 - t1half) / (t2 - t1)) * E_f1
#         + f2 * E_f2
#         + f3 * E_f3)

#     return np.exp(-1.0 * E_bar) # Return B
