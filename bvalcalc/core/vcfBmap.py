from bvalcalc.utils.load_vcf import load_vcf


def vcfBmap(args, vcf_path):    

    chrom_array, pos_array = load_vcf(vcf_path)
    print(pos_array, chrom_array)
    print("gottem in vcfBmap")


