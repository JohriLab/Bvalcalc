## Figure 2



## Figure 4
Bvalcalc --generate_params dromel_cds
Bvalcalc --download_sample_data

Bvalcalc \
--region chr2R:8000000-9000000 \
--params ./DroMel_Cds_Params.py \
--bedgff ./cds_noX.bed \
--rec_map ./dmel_comeron_recmap.csv \
--plot 1Mb_B.png

# Mean B of neutral sites across specified region: 0.8925140916007113
# Plot saved to 1Mb_B.png
# = B value calculated in 6.67 seconds. = = =
