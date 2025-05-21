from bvalcalc.utils.load_vcf import load_vcf
from bvalcalc.utils.load_Bmap import load_Bmap


def vcfBmap(args, vcf_path):    

    chrom_array, pos_array = load_vcf(vcf_path)
    print(pos_array, chrom_array)

    chromosomes, positions, b_values = load_Bmap(file_path = args.Bmap)
    print(chromosomes, positions, b_values)
    print("gottem in vcfBmap")


