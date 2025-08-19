# Sample Data for Bvalcalc

This directory contains sample data files from the literature that are compatible with Bvalcalc.

## Available Files

### Parameter Files

- `DroMel_Cds_Params.py` - Drosophila melanogaster parameters for coding sequences
- `DroMel_Utr_Params.py` - Drosophila melanogaster parameters for UTR regions
- `DroMel_Phastcons_Params.py` - Drosophila melanogaster parameters for PhastCons conserved regions

### Annotation Files

- `cds_noX.bed` - Coding sequence annotations (BED format)
- `utr_noX.bed` - UTR region annotations (BED format)
- `phastcons_noX.bed` - PhastCons conserved region annotations (BED format) (Siepel et al 2005)

### Recombination Maps

- `dmel_comeron_recmap.csv` - Drosophila melanogaster recombination map (Comeron et al 2012)

## Usage

These files can be used with Bvalcalc commands. For example:

```bash
# Download sample data
bvalcalc --download_sample_data

# Use with Bvalcalc commands
bvalcalc --genome --params DroMel_Cds_Params.py --bedgff cds_noX.bed --rec_map dmel_comeron_recmap.csv
```
