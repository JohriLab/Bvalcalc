#calculate_B_analytically_Eq3_mine_demography_multiple_elements
#This script is to get B values across a genomic element with multiple functional elements.
#Currently the output is in terms of an average of a sliding window

import sys
import math
import numpy as np
import csv

#Define variables and constants:
out_folder="droso_single_exon_gc_10kb_decline10x"
g = 0.5*1e-8 #1e-8 #rate of gene conversion
tract_len=440 #mean tract length of gene conversion in base pairs
r = 0.5*1.0*1e-8 #rate of recombination
u = 3.0*1e-9 #(*Mutation rate*)
l = 10000
U = l*u
chr_start = 1
chr_end = 2200

blockstart = []
blockend = []

with open('../exampleData/test.bed', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if len(row) >= 2:
            blockstart.append(int(row[1]))
            blockend.append(int(row[2]))

print(blockend) #todelete

blockstart = 1000 #todelete
blockend = 2000 #todelete

# blockstart = oneline

#Parameters of genome architecture:
#!!!Read in a bed file with positions of functional elements and the last position (i.e., where the genome ends)#!!!
# !!!Save in memory, the positions of the selected sites. #!!! No need to save the positions of the neutral sites (might take upp too much memory).
# !!!last_position = 10000.0 #(*Full length of the chromosome; the last position. Get this from the user or the input file.!!!*)

#Parameters of instantaneous change in demographic history:
Nanc = 1e6 #(Ancestral population size)
Ncur = 1e5 #(Current population size)
TIME="T_1" #T_0_1/T_0_5/T_1
time_of_change=1.0 #0.1/0.5/1(This is the time of change in 2Ncur generations in the past.)

#Parameters of the DFE:
DFE="DFE3"
f0 = 0.1 #(*Proportion of effectively neutral mutations with 0 <= |2Nes| < 1 *)
f1 = 0.1 #(*Proportion of weakly deleterious mutations with 1 <= |2Nes| < 10 *)
f2 = 0.1 #(*Proportion of moderately deleterious mutations with 10 <= |2Nes| < 100 *)
f3 = 0.7 #(*Proportion of strongly deleterious mutations with |2Nes| >= 100 *)
#(*Note that the number of classes can easily be increased to whatever is required to approximate the continuous DFE *)
h = 0.5 #(* dominance coefficient *)
gamma_cutoff = 2.0 #5.0
s_window_size = 100

#Constants that we do not need from the user:
pi_anc = 4*Nanc*u #(*Expected nucleotide diversity under neutrality*)
#(*Now we define the boundaries of the fixed intervals over which we will integrate. The number of bins can be a user parameter, if we like. *)
t0 = 0.0
t1 = h*(1/(2*Nanc))
t1half = h*(gamma_cutoff/(2*Nanc)) #(* This is the cut-off value of 2Nes=5. This derivation assumes that all mutations with 2Nes<5 will not contribute to BGS *)
t2 = h*(10/(2*Nanc))
t3 = h*(100/(2*Nanc))
t4 = h*1.0

neu_sites = []
for site in range(chr_start, chr_end + 1):
    if not (blockstart <= site <= blockend):
        neu_sites.append({
            "pos": site,
            "div1000": site / 1000
        })

print(neu_sites[-1])
sys.exit()

#!!!def get_distance_to_functional_element (posn, element)#!!! there is some flexibility in how to write this part.
#!!!    return()

#calculate the quantities "a" and "b" which are constants that depend on the recombination and gene conversion rate and also the distance between the focal site and the functional element.
def calculate_a_and_b(distance_to_element, length_of_element):
    C = (1.0 - math.exp(-2.0*r*distance_to_element))/2.0
    if g==0:
        a = C
        b = C + r*length_of_element
    elif g > 0:
        if distance_to_element+length_of_element < 0.5*tract_len:#The 0.5 is currently aribtrary. It's just about when the approximation of 1-exp(-x)~x holds. 
            #print ("accounting for small-distance gene conversion")
            a = (C + (g*distance_to_element))
            b = (C + r*length_of_element + (g*(distance_to_element+length_of_element)))
        else:
            #print ("accounting for large-distance gene conversion")
            a = g*tract_len + C
            b = g*tract_len + r*length_of_element + C
    return (a, b)

#calculate the exponent using previously computed values of "a" and "b"
def calculate_exponent(t_start, t_end, distance_to_element, length_of_element):
    t_tmp = calculate_a_and_b(distance_to_element, length_of_element)
    a = t_tmp[0]
    b = t_tmp[1]
    E1 = ((U*a)/((1-a)*(a-b)*(t_end-t_start))) * math.log((a+(t_end*(1-a)))/(a + (t_start*(1-a))))
    E2 = -1.0*((U*b)/((1-b)*(a-b)*(t_end-t_start)))*math.log((b + ((1-b)*t_end))/(b + ((1-b)*t_start)))
    E = E1 + E2
    return (E)

#Calculate B due to a single functional element at the focal site. Here we sum over the DFE.
def calculate_B(distance_to_element, length_of_element):
    E_f1 = calculate_exponent(t1half, t2, distance_to_element, length_of_element) 
    E_f2 = calculate_exponent(t2, t3, distance_to_element, length_of_element)
    E_f3 = calculate_exponent(t3, t4, distance_to_element, length_of_element)
    E_bar = f0*0.0 + f1*((t1half-t1)/(t2-t1))*0.0 + f1*((t2-t1half)/(t2-t1))*E_f1 + f2*E_f2 + f3*E_f3
    B = math.exp(-1.0*E_bar)
    return (B)

#calculate an average B value over the window with coordinates win_start - win_end
def calculate_Banc_window(win_start, win_end):
    i=win_start
    B_sum = 0.0
    s_tot = 0
    while i <= win_end:
        #!!!if position i is a neutral site: #!!! Check if this site is neutral. Calculate B if it is neutral. Otherwise do not compute B.
            #calculating B at site i from all genomic elements:
        #!!!    B_i = 1.0
        #!!!    For each functional element in the genome: #!!!
        #!!!        B_i = B_i * float(calculate_B(distance_to_element[i], length_of_element[i])) #!!! give the nearest distance from position i to the given element and the length of that element
       #!!!     B_sum = B_sum + B_i
        #!!!    s_tot = s_tot + 1
        i = i + 1
    if s_tot > 0: 
        return(B_sum/s_tot)
    else: #If all sites in this window are selected, then return NA.
        return ("NA")

#Gets the value of B in the current population as a function of B calcualted in the ancestral population (assumed to be in equilibrium).
#As in, this is where we account for a simple single-size change in N
def get_Bcur(Banc):
    R = float(Nanc)/float(Ncur)
    Bcur = (Banc*(1.0 + (R-1.0)*math.exp((-1.0*time_of_change)/Banc))) / (1 + (R-1.0)*math.exp(-1.0*time_of_change))
    return (Bcur)

                                                                                
#Getting an average pi over a window:
result = open("/work/users/p/j/pjohri/BGSCalculator/mytheory/" + out_folder + "/Bvalues_" + DFE + "_" + TIME + "_windowsize_" + str(s_window_size) + ".txt", 'w+')
result.write("start" + '\t' + "end" + '\t' + "B" + '\n')
posn_start = 1
while posn_start <= chr_end:
    result.write(str(posn_start) + '\t' + str(posn_start+s_window_size) + '\t' + str(get_Bcur(calculate_Banc_window(posn_start, posn_start+s_window_size))) + '\n')
    posn_start = posn_start + s_window_size
result.close()

print("done")


