import numpy as np
from core.calculateB import calculateB_unlinked
import sys

def calc_B_from_other_chromosomes(unlinked_L):
    print("Gacha")
    unlinked_B = calculateB_unlinked(unlinked_L)
    print(unlinked_B)
    # sys.exit()