from core.chromBcalc import chromBcalc

def genomeBcalc(args):    
    output_data, block_ranges = chromBcalc(args)
    return output_data, block_ranges