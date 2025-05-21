from bvalcalc.utils.load_vcf import load_vcf
from bvalcalc.utils.load_Bmap import load_Bmap
import numpy as np


def vcfBmap(args, vcf_path):    

    vcf_chromosomes, vcf_positions = load_vcf(vcf_path)
    print(vcf_chromosomes, vcf_positions)

    bmap_chromosomes, bmap_positions, b_values = load_Bmap(file_path = args.Bmap)
    print(bmap_chromosomes, bmap_positions, b_values)

    max_bmap_pos = bmap_positions[-1]
    above_count  = np.sum(vcf_positions > max_bmap_pos)
    if above_count > 0:
        print(f"WARNING: {above_count} VCF positions are above the max B-map position ({max_bmap_pos})")


    print("gottem in vcfBmap")


