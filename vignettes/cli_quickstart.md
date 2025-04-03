# 🔧 Bvalcalc CLI Quickstart

Welcome to **Bvalcalc**, your command-line companion for calculating B-values! This quickstart guide will show you how to get up and running in just a few steps.

## 📦 Installation

You can install Bvalcalc via pip:

```bash
pip install bvalcalc
```

Okay now lets see what a section of the Drosophila genome might look like

```bash
./Bvalcalc.py --genome --pop_params drosophila_params.py --bedgff_path exampleData/dmel6_2R_genes.csv --plot_output --calc_start 12000000 --calc_end 13000000
```
