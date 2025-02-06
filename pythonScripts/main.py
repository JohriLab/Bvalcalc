#!/usr/bin/env python3

# ./Bcalc_multiElements.py --chr_start 1 --chr_end 25254535 --file_path ../exampleData/dmel6_2R_genes.csv

#calculate_B_analytically_Eq3_mine_demography_multiple_elements
#This script is to get B values across a genomic element with multiple functional elements.
#Currently the output is in terms of an average of a sliding window
import sys
from helperScripts.parseArgs import parseArgs
from helperScripts.genomeBcalc import genomeBcalc
from helperScripts.regionBcalc import regionBcalc
from helperScripts.plotBasic import plotBasic
from helperScripts.calculateB import calculateB
import time

#CLI handling
def main():
    start_time = time.time()

    args = parseArgs() 
    # genomeBcalc(args)
    regionBcalc(args)
    # plotBasic(args)
    # print("4k test", calculateB(4000, 10000))
    # print("1.25k test", calculateB(1250, 10000))
    

    print(f"Script completed in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()

sys.exit()