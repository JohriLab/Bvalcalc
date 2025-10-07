## All results are using Bvalcalc v1.0.0
## We have these results as unit and end to end tests so compatibility should be maintained going forward

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

## Table SA.1
# Note that this was done using the Python REPL to access the API functions
from Bvalcalc import get_params, calculateB_unlinked, calculateB_linear, calculateB_hri
params = get_params("tests/testparams/InterfParams.py", constant_dfe = True)
calculateB_linear(distance_to_element = 0, length_of_element = 10000, params = params)
calculateB_hri(distant_B = 1, interfering_L = 10000, params = params)
##

## Table SA.2
# Note that this was done using the Python REPL to access the API functions
from Bvalcalc import get_params, calculateB_unlinked, calculateB_linear, calculateB_hri
params = get_params("tests/testparams/nogcBasicParams.py")
calculateB_linear(distance_to_element = 500, length_of_element = 10000, params = params)
calculateB_hri(distant_B = 1, interfering_L = 10000, params = params)
##