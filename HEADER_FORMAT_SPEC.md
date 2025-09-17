# Bvalcalc Header Format Specification

## Overview

This specification defines how Bvalcalc handles header comments in input and output files. Headers are optional and use `#` comments to provide metadata, warnings, and input commands.

## Header Format

### Basic Structure

```
# Bvalcalc Header Format v1.0
# Generated: 2024-01-15 14:30:25
# Command: bvalcalc --genome --params DroMel_Cds_Params.py --bedgff cds_noX.bed --rec_map dmel_comeron_recmap.csv --out output.csv
#
# WARNING: Recombination map doesn't cover centromeric regions
# WARNING: Using default rate (1.0) for missing regions
#
# Parameters:
#   - r: 1.2e-8
#   - g: 1.0e-6
#   - Nanc: 1000000
#
# Input Files:
#   - BED/GFF: cds_noX.bed (Drosophila CDS annotations)
#   - Rec Map: dmel_comeron_recmap.csv (Comeron et al. 2012)
#   - Params: DroMel_Cds_Params.py
#
# Data: Chromosome,Start,B
chr2L,1,0.987654
chr2L,1000,0.985432
...
```

### Header Types

#### 1. **File Identification**

```
# Bvalcalc Header Format v1.0
# Generated: 2024-01-15 14:30:25
# File Type: B-map
```

#### 2. **Command Information**

```
# Command: bvalcalc --genome --params DroMel_Cds_Params.py --bedgff cds_noX.bed --rec_map dmel_comeron_recmap.csv --out output.csv
# Mode: genome
# Parameters File: DroMel_Cds_Params.py
```

#### 3. **Warnings**

```
# WARNING: Recombination map doesn't cover centromeric regions
# WARNING: Using default rate (1.0) for missing regions
# WARNING: No chromosome size found for 'chrM' in chr_sizes.csv
```

#### 4. **Input Files**

```
# Input Files:
#   - BED/GFF: cds_noX.bed (Drosophila CDS annotations)
#   - Rec Map: dmel_comeron_recmap.csv (Comeron et al. 2012)
#   - Params: DroMel_Cds_Params.py
#   - Chr Sizes: dmel6_chr_sizes.csv
```

#### 5. **Parameters**

```
# Parameters:
#   - r: 1.2e-8 (recombination rate)
#   - g: 1.0e-6 (gene conversion rate)
#   - Nanc: 1000000 (ancestral population size)
#   - h: 0.5 (dominance coefficient)
```

#### 6. **Data Description**

```
# Data: Chromosome,Start,B
# Format: CSV with chromosome name, start position, B-value
```

## Implementation Details

### Header Parsing

- Headers start with `#` and continue until the first non-comment line
- Headers are preserved when reading files
- Headers are automatically generated when writing files
- Existing files without headers continue to work (backward compatibility)

### Header Generation

- **Input files**: Headers are preserved and passed through
- **Output files**: Headers are automatically generated with:
  - File identification
  - Command used
  - Warnings encountered
  - Input files used
  - Parameters used
  - Data format description

### File Types Supported

1. **BED/GFF files** - Annotation files
2. **CSV files** - Recombination maps, gene conversion maps, chromosome sizes
3. **B-map files** - Output files with B-values
4. **VCF files** - Variant files (if applicable)

## Examples

### Input File with Header

```
# Bvalcalc Header Format v1.0
# Generated: 2024-01-15 14:30:25
# File Type: BED
# Description: Drosophila CDS annotations (excluding X chromosome)
# Source: FlyBase release 6.32
#
# Data: Chromosome,Start,End
chr2L,1000,2000
chr2L,3000,4000
...
```

### Output File with Header

```
# Bvalcalc Header Format v1.0
# Generated: 2024-01-15 14:30:25
# Command: bvalcalc --genome --params DroMel_Cds_Params.py --bedgff cds_noX.bed --rec_map dmel_comeron_recmap.csv --out output.csv
# Mode: genome
#
# WARNING: Recombination map doesn't cover centromeric regions
# WARNING: Using default rate (1.0) for missing regions
#
# Input Files:
#   - BED/GFF: cds_noX.bed (Drosophila CDS annotations)
#   - Rec Map: dmel_comeron_recmap.csv (Comeron et al. 2012)
#   - Params: DroMel_Cds_Params.py
#
# Parameters:
#   - r: 1.2e-8
#   - g: 1.0e-6
#   - Nanc: 1000000
#
# Data: Chromosome,Start,B
chr2L,1,0.987654
chr2L,1000,0.985432
...
```

## Backward Compatibility

- Files without headers continue to work unchanged
- Headers are optional - existing functionality is preserved
- New header format is additive, not replacing existing functionality
