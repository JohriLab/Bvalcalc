import numpy as np

def load_vcf(vcf_path):
    print(f"loading VCF/CSV: {vcf_path}")
    chromosomes = []
    positions   = []

    with open(vcf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Detect comma-separated two-column CSV vs. whitespace-separated VCF
            if ',' in line:
                parts = [p.strip() for p in line.split(',')]
                # only accept exactly two columns
                if len(parts) != 2:
                    continue
            else:
                parts = line.split()
                # need at least chromosome and POS
                if len(parts) < 2:
                    continue

            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                # bad POS field, skip
                continue

            chromosomes.append(chrom)
            positions.append(pos)

    chrom_array = np.array(chromosomes, dtype='<U20')
    pos_array   = np.array(positions,   dtype=np.int64)
    return chrom_array, pos_array
