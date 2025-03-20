import numpy as np
from ExampleParams import Nanc, Ncur, time_of_change

# #Gets the value of B in the current population as a function of B calcualted in the ancestral population (assumed to be in equilibrium).
# #As in, this is where we account for a simple single-size change in N

def get_Bcur(Banc):
    R = float(Nanc)/float(Ncur)
    Bcur = (Banc*(1.0 + (R-1.0)*np.exp((-1.0*time_of_change)/Banc))) / (1 + (R-1.0)*np.exp(-1.0*time_of_change)) ##Currently time_of_change is hardcoded, need to change to variable
    # Bnew = 1 - (1-Banc)*(Ncur/Nanc)
    # Bcur = Bnew + (Banc - Bnew) * np.exp(-time_of_change/(2*Ncur))
    return (Bcur)
