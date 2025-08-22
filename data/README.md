# Drosophila CDS Sample Data for Bvalcalc

This directory contains Drosophila melanogaster CDS sample data files from the literature that are compatible with Bvalcalc.

## Available Files

### Annotation Files

- `cds_noX.bed` - Coding sequence annotations (BED format)

### Recombination Maps

- `dmel_comeron_recmap.csv` - Drosophila melanogaster recombination map (Comeron et al 2012)

## Usage

These files can be used with Bvalcalc commands. For example:

```bash
# Download sample data to current directory
bvalcalc --download_sample_data

# Download sample data to specific directory
bvalcalc --download_sample_data --dir /path/to/directory

# Use with Bvalcalc commands
bvalcalc --genome --params your_params.py --bedgff cds_noX.bed --rec_map dmel_comeron_recmap.csv
```
