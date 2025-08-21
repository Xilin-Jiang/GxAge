##### ! As all these analysis request individual level data which cannot be shared,
##### ! we shared our summary data and plotting.
##### ! We provide a simulated example PRS + disease traits so that people can replicate
##### ! the analysis by formatting their data in the same format.
##### ! these simulated data are created by permuting and add noise to real data so they have realistic distribution.

##########################################################
# example of estimating age dependent prediction accuracy
##########################################################

## formatting your data:
# indiviudal_id
# disease: 1 if disease is observed, 0 if not
# age_diag: age point of the disease diagnosis; will be ignored if for rows with disease == 0
# censor_age: the last observation of the individual
# RiskScore: the risk score of interest
# age_risk_measure: age when the risk score was estimated, will be regressed out from risk score; if not available, just put 0 for all
# genetic_sex: genetic sex of the individual, will be regressed out from risk score; if not available, just put 0 for all
# PC1:PC5: additional covariates that should be regressed out from risk score; if not available, just put 0 for all
example_data <- fread("Summary_data/example_RiskScore_age_data.txt")

source("GxAge_functions.R")
results <- Prediction_R2_Slope_per10years(example_data, flag_normal_transform_RiskScore = T)
print(paste0("per 10 year change: ", results[[2]]$incident_R2_1to1_slope_10year,
             " z-score = ", results[[2]]$slope_zscore))

