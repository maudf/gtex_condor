#!/bin/bash

for i in {1..5}
do
    Rscript resample_GO_results_shared.R lung $i 200 >log_resampling/log_LNG_$i.txt &
    Rscript resample_GO_results_shared.R whole_blood $i 200 >log_resampling/log_WBL_$i.txt &
done

for i in {1..10}
do
    Rscript resample_GO_results_shared.R artery_aorta $i 100 >log_resampling/log_ATA_$i.txt &
    Rscript resample_GO_results_shared.R artery_tibial $i 100 >log_resampling/log_ATT_$i.txt &
    Rscript resample_GO_results_shared.R esophagus_muscularis $i 100 >log_resampling/log_EMS_$i.txt &
    Rscript resample_GO_results_shared.R heart_left_ventricle $i 100 >log_resampling/log_HRV_$i.txt &
done

for i in {1..20}
do
    Rscript code/resample_GO_results_shared.R adipose_subcutaneous $i 50 >code/log_resampling/log_ADS_$i.txt &
    Rscript code/resample_GO_results_shared.R cells_transformed_fibroblasts $i 50 >code/log_resampling/log_FIB_$i.txt &
    Rscript code/resample_GO_results_shared.R esophagus_mucosa $i 50 >code/log_resampling/log_EMC_$i.txt &
    Rscript code/resample_GO_results_shared.R muscle_skeletal $i 50 >code/log_resampling/log_SKM_$i.txt &
    Rscript code/resample_GO_results_shared.R nerve_tibial $i 50 >code/log_resampling/log_TNV_$i.txt &
    Rscript code/resample_GO_results_shared.R skin $i 50 >code/log_resampling/log_SKN_$i.txt &
    Rscript code/resample_GO_results_shared.R thyroid $i 50 >code/log_resampling/log_THY_$i.txt &
done
