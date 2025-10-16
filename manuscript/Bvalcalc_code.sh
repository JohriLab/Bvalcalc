## All results are using Bvalcalc v1.0.0 CLI and Python API
## We have these results as unit and end to end tests so compatibility should be maintained going forward
## Executed from root directory of Bvalcalc repository

## Figure 2
# 2a Basic model
Bvalcalc --gene --params tests/testparams/nogcBasicParams.py

# 2b Gene conversion
Bvalcalc --gene --params tests/testparams/gcBasicParams.py

# 2c Population expansion
Bvalcalc --gene --params tests/testparams/ExpandParams_5N_1T.py --pop_change

# 2d Selfing
Bvalcalc --gene --params tests/testparams/SelfParams_0.9S_0.5h.py
##

## Supplementary Figures
Bvalcalc --gene --params tests/testparams/ExpandParams_5N_0.2T.py --pop_change
Bvalcalc --gene --params tests/testparams/ContractParams_5N_1T.py --pop_change
Bvalcalc --gene --params tests/testparams/ContractParams_5N_0.2T.py --pop_change
# See Figure3.py for 200kb plot of basic model
##

## Figure 4
Bvalcalc --generate_params dromel_cds
Bvalcalc --download_sample_data

Bvalcalc \
--region chr2R:8000000-9000000 \
--params ./DroMel_Cds_Params.py \
--bedgff ./cds_noX.bed \
--rec_map ./dmel_comeron_recmap.csv \
--plot 1Mb_B.png
##

## Table 1
# Note that I dynamically changed the UnlinkedParams.py DFE and re-ran the following command for each condition
# This could've been done more efficiently using the calculateB_unlinked function in the Python API
Bvalcalc --region chr_neutral:1-1 --params tests/testparams/UnlinkedParams.py --bedgff tests/testfiles/200kb_unlinked.csv
Bvalcalc --region chr_neutral:1-1 --params tests/testparams/UnlinkedParams.py --bedgff tests/testfiles/200kb_unlinked.csv --constant_dfe
##

## Table SA.1 B' with constant DFE
# This was done using the Python REPL to access the API functions
# Note that I dynamically changed the HriParams.py DFE and re-ran the following command for each condition
# Currently need to flush the cache each time to get the correct results
from Bvalcalc import get_params, calculateB_unlinked, calculateB_linear, calculateB_hri
params = get_params("tests/testparams/HriParams.py", constant_dfe = True)
calculateB_linear(distance_to_element = 0, length_of_element = 10000, params = params)
calculateB_hri(distant_B = 1, interfering_L = 10000, params = params)
##

## Table SA.2 B' with variable DFE
# This was done using the Python REPL to access the API functions
# Note that I dynamically changed the HriParams.py DFE and re-ran the following command for each condition
# Currently need to flush the cache each time to get the correct results
from Bvalcalc import get_params, calculateB_unlinked, calculateB_linear, calculateB_hri
params = get_params("tests/testparams/HriParams.py")
calculateB_linear(distance_to_element = 0, length_of_element = 10000, params = params)
calculateB_hri(distant_B = 1, interfering_L = 10000, params = params)
##

## Bmap code
# See the manuscript for how the bedgffs and recmaps were sourced
# HomSap
Bvalcalc --genome --params ./HomSap_Cds_Params.py --bedgff ../bedgffs/HomSap_cds.bed --rec_map ../recmaps/HomSap_recmap.csv --gc_map ../recmaps/HomSap_gcmap.csv --pop_change --out ./homsap_cds_bmap.csv --out_binsize 1000 --verbose  --chr_sizes ./chromsizes.txt
Bvalcalc --genome --params ./HomSap_Phastcons_Params.py --bedgff ../bedgffs/HomSap_phastcons.bed --rec_map ../recmaps/HomSap_recmap.csv --gc_map ../recmaps/HomSap_gcmap.csv --pop_change --prior_Bmap ./homsap_cds_bmap.csv --out ./homsap_cdsphastcons_bmap.csv --out_binsize 1000 --verbose --chr_sizes ./chromsizes.txt
#
# DroMel
Bvalcalc --genome --gamma_dfe \
    --pop_params DroMel_Cds_Params.py \
    --bedgff_path ../bedgffs/DroMel_cds.bed \
    --out dromel_cds_bmap.csv \
    --out_binsize 1000 --rec_map ../recmaps/DroMel_recmap.csv
Bvalcalc --genome --gamma_dfe \
    --pop_params DroMel_Phastcons_Params.py \
    --bedgff_path ../bedgffs/DroMel_phastcons.bed \
    --prior_Bmap dromel_cds_bmap.csv \
    --out bmaps/highmutrate/phastcons_CDS_B.csv \
    --out_binsize 1000 --rec_map ../recmaps/DroMel_recmap.csv
Bvalcalc --genome --gamma_dfe \
    --pop_params DroMel_Utr_Params.py \
    --bedgff_path ../bedgffs/DroMel_utr.bed \
    --prior_Bmap dromel_cds_phastcons_bmap.csv \
    --out dromel_cds_phastcons_utr_bmap.csv \
    --out_binsize 1000 --rec_map ../recmaps/DroMel_recmap.csv
#
# AraTha
Bvalcalc --genome --params ./AraTha_Cds_Params.py --bedgff ../bedgffs/AraTha_cds.gff3 --rec_map ../recmaps/AraTha_recmap.csv --pop_change --out ./aratha_cds_bmap.csv --out_binsize 1000 --verbose
#