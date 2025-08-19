# Sample Data for Bvalcalc

This directory contains sample data files for testing and learning with Bvalcalc.

## Available Files

### Parameter Files

- `DroMel_Cds_Params.py` - Drosophila melanogaster parameters for coding sequences
- `DroMel_Utr_Params.py` - Drosophila melanogaster parameters for UTR regions
- `DroMel_Phastcons_Params.py` - Drosophila melanogaster parameters for PhastCons conserved regions

### Annotation Files

- `cds_noX.bed` - Coding sequence annotations (BED format)
- `utr_noX.bed` - UTR region annotations (BED format)
- `phastcons_noX.bed` - PhastCons conserved region annotations (BED format)

### Recombination Maps

- `dmel_comeron_recmap.csv` - Drosophila melanogaster recombination map

## Usage

These files can be used with Bvalcalc commands. For example:

```bash
# Download sample data
bvalcalc --download_sample_data

# Use with Bvalcalc commands
bvalcalc --genome --params DroMel_Cds_Params.py --bedgff cds_noX.bed --rec_map dmel_comeron_recmap.csv
```
