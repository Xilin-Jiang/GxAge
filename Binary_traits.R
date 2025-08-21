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

##########################################################
# Age-dependent prediction accuracy of PRS
##########################################################
# the methods are as above, we can't share individual-level data
PRS_binary_cox_effect_each_age <- fread("Summary_data/SuppFigure4.csv")

curve_cox_effect <- ggplot(PRS_binary_cox_effect_each_age, mapping = aes(x = age_bin_idx, y = cox_effects, group = binary_name)) +
  geom_line(aes(group = binary_name), linewidth = 0.5, color =blue, linetype = 2) +
  geom_pointrange(aes(ymin=cox_effects - cox_se, ymax=cox_effects+ cox_se), size = 0.5,  color = blue, alpha=0.7) +
  scale_fill_manual(values = rep(blue, 5)) +
  # coord_flip() +
  theme_bw(base_size = 20)+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + labs(x = "", y = expression(paste(log~"(Hazard Ratio)"))) +
  scale_y_continuous(labels = scales::percent, ,limits = c(0,NA)) +
  facet_grid( ~ binary_name)
ggsave(paste0("Figures/curve_cox_effect_acrossAgeBins.pdf"), curve_cox_effect, width = 12, height = 4)

# prevalent vs. incident
slope_results <- fread(sep = "\t", paste0("Summary_data/PrevelantAssociation_PRS_slope_and_popLevelEstiamtes_20250403.txt"))

slope_R2 <- slope_results %>%
  mutate(R2_bonforonni = (1 - pnorm(abs(observedR2diseasResidual_slope_per_year/bootstrapping_se_observedR2diseasResidual_slope_per_year)) ) * 2  * length(binary_names_list$V1)) %>%
  mutate(prevelant_R2_percent_per_10year = observedR2diseasResidual_slope_per_year * 10 / all_observed_R2_diseaseResidual,
         prevelant_R2_percent_per_10year_se = bootstrapping_se_observedR2diseasResidual_slope_per_year * 10 / all_observed_R2_diseaseResidual, R2_significant = R2_bonforonni < 0.1) %>%
  select(binary_name, prevelant_R2_percent_per_10year, prevelant_R2_percent_per_10year_se,
         R2_bonforonni, all_observed_R2, R2_significant) %>%
  arrange(prevelant_R2_percent_per_10year)

# load incident case prediction accuracy
incident_slope_results <- fread(sep = "\t", paste0("Summary_data/slope_and_popLevelEstiamtes_CoxLinearLogisticR2Estimates_20250310.txt"))

incident_slope_R2 <- incident_slope_results %>%
  mutate(R2_bonforonni = (1 - pnorm(abs(R2_slope_per_year/bootstrapping_se_R2_slope_per_year)) ) * 2  * length(binary_names_list$V1)) %>%
  mutate(R2_percent_per_10year = R2_slope_per_year * 10 / all_observed_R2,
         R2_percent_per_10year_se = bootstrapping_se_R2_slope_per_year * 10 / all_observed_R2, R2_significant = R2_bonforonni < 0.1) %>%
  select(binary_name, R2_percent_per_10year, R2_percent_per_10year_se,
         R2_bonforonni, all_observed_R2, R2_significant) %>%
  setNames(paste0('incident_', names(.)))

df_incidentVSprevelant_change <- slope_R2 %>%
  left_join(incident_slope_R2, by = c("binary_name" = "incident_binary_name"))
plt_incidentVSprevelant_change <- ggplot(data=df_incidentVSprevelant_change) +
  geom_pointrange(aes(
    y=incident_R2_percent_per_10year, ymin=incident_R2_percent_per_10year - incident_R2_percent_per_10year_se,
    ymax = incident_R2_percent_per_10year + incident_R2_percent_per_10year_se,  x=prevelant_R2_percent_per_10year), color = red,
    size = 1, alpha = 1) +
  geom_pointrange(aes(
    y=incident_R2_percent_per_10year, x=prevelant_R2_percent_per_10year,
    xmin = prevelant_R2_percent_per_10year - prevelant_R2_percent_per_10year_se,
    xmax = prevelant_R2_percent_per_10year + prevelant_R2_percent_per_10year_se), color = red, shape = 21,
    size = 1, alpha = 1) +
  scale_fill_gradient(limits = c(0,0.7), low = "white", high = red) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(y = expression(paste( ~Delta, R^2, "(PRS, incident cases) w.r.t age (per 10 years)")), x = expression(paste( ~Delta, R^2, "(PRS, prevelant cases) w.r.t age (per 10 years)"))) +
  geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.9,0.3))+
  scale_x_continuous(labels = scales::percent, limits = c(-0.9,0.3)) +
  ggrepel::geom_label_repel(data = filter(df_incidentVSprevelant_change),
                            aes(y=incident_R2_percent_per_10year, x=prevelant_R2_percent_per_10year, label = binary_name)
  )

ggsave(paste0("Figures/PRS_incidentVSprevelant_20250404.pdf"), plt_incidentVSprevelant_change, width = 6, height = 6)


#############################################
# EA liability-threshold model vs. PRS R2
#############################################
pop_h2g_predictedliability_results <- fread( paste0("Summary_data/predLiab_h2g_results_all_20250311.txt"))
h2g_predictliability_change_test <- fread(file = "Summary_data/PredictedLiability_h2g_change_with_age_qtrait_likelihoodRatio_test_20250122.txt", sep = "\t")
plot_h2g_change <- h2g_predictliability_change_test %>%
  left_join(pop_h2g_predictedliability_results, by = c("trait_list" = "quant_name")) %>%
  mutate(h2g_pop_normalised_percent_per_10year = trait_slope_per_year * 10 / h2g_mean,
         h2g_pop_normalised_percent_per_10year_se = slope_se * 10 / h2g_mean, significant = bonforonni < 0.1) %>%
  arrange(h2g_pop_normalised_percent_per_10year) %>%
  select(trait_list, h2g_pop_normalised_percent_per_10year, h2g_pop_normalised_percent_per_10year_se,
         h2g_mean)

slope_results <- fread(sep = "\t", paste0("Summary_data/QRSmatchPRSaccuracy_binary_R2_slope_and_popLevelEstiamtes_20250415.txt"))
slope_R2 <- slope_results %>%
  mutate(binary_name = paste0(binary_name, "_liability"),
         QRS_R2_percent_per_10year = R2_slope_per_year * 10 / all_observed_R2,
         QRS_R2_percent_per_10year_se = mc_se_R2_slope_per_year * 10 / all_observed_R2,
         PRS_R2_percent_per_10year = PRS_R2_slope_per_year * 10 / all_PRS_observed_R2,
         PRS_R2_percent_per_10year_se = mc_se_PRS_R2_slope_per_year * 10 / all_PRS_observed_R2) %>%
  select(binary_name, QRS_R2_percent_per_10year, QRS_R2_percent_per_10year_se,
         PRS_R2_percent_per_10year, PRS_R2_percent_per_10year_se, all_observed_R2,all_PRS_observed_R2) %>%
  arrange(QRS_R2_percent_per_10year)

binary_names_list <- fread("Summary_data/binarytraits_list.txt", header = F) %>%
  mutate(liability_name = paste0(V1, "_liability"))


###### ! in this section we assume sampling correlation between h2g and QRS prediction effects is 0 as bootstraping population to get h2g + QRS estimates are
# computationally infeasible. H2g and QRS prediction effects are approximately independent during the sampling process, as
#  QRS and h2g are estimated using different data sets (QRS uses diseases, h2g uses continuous measurements + genetics)
###### ! in all other prediction analysis, we perform monte-carlo sampling of individuals to get point estimates and se.

QRS_h2g_liability_each_age <- list()
for(idx in 1:length(binary_names_list$liability_name)){
  liability_name <- binary_names_list$liability_name[idx]
  binary_name <- binary_names_list$V1[idx]
  print(liability_name)

  QRS_R2_percent_per_10year <- slope_R2 %>%
    filter(binary_name == !!liability_name) %>%
    pull(QRS_R2_percent_per_10year)
  QRS_R2_percent_per_10year_se <- slope_R2 %>%
    filter(binary_name == !!liability_name) %>%
    pull(QRS_R2_percent_per_10year_se)

  h2g_pop_normalised_percent_per_10year <- plot_h2g_change %>%
    filter(trait_list == !!liability_name) %>%
    pull(h2g_pop_normalised_percent_per_10year)
  h2g_pop_normalised_percent_per_10year_se <- plot_h2g_change %>%
    filter(trait_list == !!liability_name) %>%
    pull(h2g_pop_normalised_percent_per_10year_se)

  QRS_point_estimates_expected <- (1+QRS_R2_percent_per_10year) * (1+h2g_pop_normalised_percent_per_10year)  - 1
  # independence: QRS and h2g are estimated using different data sets (QRS uses diseases, h2g uses continuous measurements + genetics)
  mc_QRS_samples <- rnorm(100, mean = QRS_R2_percent_per_10year, sd = QRS_R2_percent_per_10year_se)
  mc_h2g_samples <- rnorm(100, mean = h2g_pop_normalised_percent_per_10year, sd = h2g_pop_normalised_percent_per_10year_se)

  QRS_expected_se <- sd( (1+mc_QRS_samples) * (1+mc_h2g_samples)  - 1 )

  QRS_h2g_liability_each_age[[idx]] <- data.table(binary_name, liability_name, QRS_point_estimates_expected, QRS_expected_se)

}
QRS_h2g_liability_each_age <- bind_rows(QRS_h2g_liability_each_age)

# plot slope
slope_GplusE <- fread(sep = "\t", paste0("Summary_data/QRSmatchPRSaccuracy_binary_R2_slope_and_popLevelEstiamtes_20250415.txt"))
plt_predictliabiliy <- slope_GplusE %>%
  mutate(PRS_R2_percent_per_10year = PRS_R2_slope_per_year * 10 / all_PRS_observed_R2,
         PRS_R2_percent_per_10year_se = mc_se_PRS_R2_slope_per_year * 10 / all_PRS_observed_R2) %>%
  select(binary_name, PRS_R2_percent_per_10year, PRS_R2_percent_per_10year_se) %>%
  left_join(QRS_h2g_liability_each_age, by = "binary_name")
plt_predictliabiliy <- plt_predictliabiliy %>%
  mutate(abs_diff_zscore = abs(PRS_R2_percent_per_10year - QRS_point_estimates_expected)/
           sqrt(PRS_R2_percent_per_10year_se^2 + QRS_expected_se^2)) %>%
  mutate(abs_diff_zscore = pmin(abs_diff_zscore, 3)) # cap it at 3

plt_QRSandh2g_vs_PRS_matched <- ggplot(data=plt_predictliabiliy) +
  geom_pointrange(aes(
    y=PRS_R2_percent_per_10year, ymin=PRS_R2_percent_per_10year - PRS_R2_percent_per_10year_se,
    ymax = PRS_R2_percent_per_10year + PRS_R2_percent_per_10year_se,  x=QRS_point_estimates_expected,
    fill = abs_diff_zscore), shape = 21,
    size = 1) +
  geom_pointrange(aes(
    y=PRS_R2_percent_per_10year, x=QRS_point_estimates_expected,
    xmin = QRS_point_estimates_expected - QRS_expected_se,
    xmax = QRS_point_estimates_expected + QRS_expected_se,
    fill = abs_diff_zscore),  shape = 21,
    size = 1) +
  scale_fill_gradient(limits = c(0,3), low = "white", high = red) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 19),
        axis.text = element_text(size = 19),
        axis.title= element_text(size = 19) ) +
  labs(y = expression(paste( ~Delta, R^2, " of PRS w.r.t age (per 10 years)")), x = expression(paste( ~Delta, R^2, " of QRS+EA w.r.t age (per 10 years)"))) +
  geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(-1.05,0.4))+
  scale_x_continuous(labels = scales::percent, limits = c(-1.05,0.4)) +
  ggrepel::geom_label_repel(data = filter(plt_predictliabiliy),
                            aes(y=PRS_R2_percent_per_10year, x=QRS_point_estimates_expected, label = binary_name),size=5
  )

ggsave(paste0("Figures/matchingR2QRS_h2gEAPA_vs_PRS_20250422.pdf"), plt_QRSandh2g_vs_PRS_matched, width = 6, height = 6)

#############################################
# EA liability-threshold model vs. PRS R2
#############################################
prediction_different_models <- fread(paste0("Summary_data/Figure7AB.csv"))

models_list <- c("average_4and5", "average_top3", "top3prediction", "liability_transformed")
models_df <- prediction_different_models %>%
  filter(predictor %in% models_list) %>%
  mutate(predictor = factor(predictor, levels = models_list))

plot_summarise <- models_df %>%
  filter(liability_R2 > 0.05) %>%
  group_by(predictor) %>%
  summarise(avg_10_year_change = mean(per_10year_change),
            avg_10_year_change_se = sqrt(sum(per_10year_change_se^2))/9,
            avg_liab_R2 = mean(liability_R2))
models_df %>%
  filter(liability_R2 < 0.05)

plt_avg_prediction_age_dependency <- ggplot(plot_summarise, mapping = aes(x = predictor, y = avg_10_year_change, fill = predictor)) +
  geom_col( position=position_dodge(0.8), width=.7, alpha = 0.8) +
  geom_errorbar(aes(ymin=avg_10_year_change - avg_10_year_change_se, ymax=avg_10_year_change+ avg_10_year_change_se), width=.2, position=position_dodge(0.8)) +
  scale_fill_manual(values = c(red, green, blue, grey)) +
  coord_flip() +
  theme_bw(base_size = 20)+
  theme(
    legend.position="none",
    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title= element_text(size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) + labs(x = "", y = expression(paste(~Delta,R^2," w.r.t. age (per 10 year)"))) +
  scale_y_continuous(labels = scales::percent, limits = c(NA, NA))
ggsave(paste0("Figures/avg_prediction_age_dependency.pdf"), plt_avg_prediction_age_dependency, width = 4, height = 4)

df_points <- prediction_different_models %>%
  filter(liability_R2 > 0.05) %>%
  mutate(percentile_liability_R2 = ntile(liability_R2, n = 10)) %>%
  mutate(percentile_liability_R2 = factor(percentile_liability_R2)) %>%
  group_by(percentile_liability_R2) %>%
  summarise(mean_10_year = mean(per_10year_change),
            mean_liability_R2 = mean(liability_R2),
            se_10_year = sd(per_10year_change)/sqrt(n()))

plt_liabR2_vs_decrease <- ggplot(df_points, mapping = aes(x = mean_liability_R2, y = mean_10_year)) +
  geom_pointrange(aes(ymin=mean_10_year - se_10_year, ymax=mean_10_year+ se_10_year) ,
                  size =0.5, color = red, alpha = 0.8) +
  # coord_flip() +
  theme(
    legend.position="none",
    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title= element_text(size = 15)
  ) + labs(x = expression(paste("Liability-scale ",R^2)), y = expression(paste(~Delta,R^2," w.r.t. age (per 10 year)"))) +
  scale_x_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_y_continuous(labels = scales::percent, limits = c(NA, NA))+
  # geom_smooth( formula = y ~ I(x-1) + 0, color = red, fill = red, method = "lm", alpha = 0.3,
  #              fullrange = T)
  geom_smooth( formula = y ~ x, color = red, fill = red, method = "lm", alpha = 0.3,
               fullrange = T) +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave(paste0("Figures/liabR2_vs_R2decrease.pdf"), plt_liabR2_vs_decrease, width = 4, height = 4)

## Above results can be obtained by applying the helper function in the example data section. However
## to increase efficiency, users can run all predictors for one bootstrapping sample, which significantly increase the speed.
## I copy below a function as an example, which can be run fo a list of predictors
compute_incident_1to1_manyPredictor <- function(PRS_liatility_df,predictor_names){
  # prediction results
  prediction_R2 <- matrix(NA, nrow = 3, ncol = length(predictor_names)) # first row is overall R2, second is prevalent R2, third is incident R2
  incident_df <- PRS_liatility_df %>%
    filter(prevelant_disease == 0)
  # sampling 1:1 case-controls
  interval_surv_data_cases <- incident_df %>%
    filter(incident_disease == 1)
  interval_surv_data_controls <- incident_df %>%
    filter(incident_disease == 0) %>%
    sample_n(size = dim(interval_surv_data_cases)[1] * 1, replace = FALSE)
  incident_1to1_data <- interval_surv_data_cases %>%
    bind_rows(interval_surv_data_controls)
  disease_resid <- resid(lm(data = incident_1to1_data, incident_disease ~ genetic_sex + PC1 + PC2 + PC3 + PC4 + PC5, na.action="na.exclude"))
  for(predict_name_idx in 1:length(predictor_names)){
    predictor <- predictor_names[predict_name_idx]
    overall_R2 <- cor(pull(PRS_liatility_df, !!predictor), PRS_liatility_df$disease, use = "pairwise.complete.obs")^2
    prev_ds <- mean(PRS_liatility_df$disease)
    overall_liab_R2 <- LiabR2(prev_ds, overall_R2)

    prevelance_R2 <- cor(pull(PRS_liatility_df, !!predictor), PRS_liatility_df$prevelant_disease, use = "pairwise.complete.obs")^2
    prev_ds <- mean(PRS_liatility_df$prevelant_disease)
    prevelance_liab_R2 <- LiabR2(prev_ds, prevelance_R2)

    # to make it comparable
    incident_joint_1to1_R2 <- cor(pull(incident_1to1_data, !!predictor), disease_resid, use = "complete")^2
    prediction_R2[,predict_name_idx] <- c(overall_liab_R2, prevelance_liab_R2, incident_joint_1to1_R2)
  }
  colnames(prediction_R2) <- predictor_names
  prediction_R2 <- data.table(metric_name = c("overall_liab_R2", "prevelance_liab_R2", "incident_joint_1to1_R2" ), prediction_R2)

  # age-bin specific incident prediction
  age_quintile <- quantile(interval_surv_data_cases$age_decimal, probs = seq(0,1,0.2))
  age_median_quintile <- quantile(interval_surv_data_cases$age_decimal, probs = seq(0.1,1,0.2))

  Age_1to2_incident_R2 <- matrix(NA, nrow = 5, ncol = length(predictor_names))
  avg_incident_age <- c()
  #se_Age_1to2_incident_R2 <- matrix(NA, nrow = 5, ncol = length(predictor_names))
  for(age_bin_idx in 1:5){
    age_start <- age_quintile[age_bin_idx]
    age_end <- age_quintile[age_bin_idx + 1]
    interval_data <-  incident_df %>%
      filter(age_decimal >= age_start, age_decimal < age_end)

    incident_age <- interval_data %>%
      filter(incident_disease == 1) %>%
      mutate(incident_age = disease_age - age_decimal)
    avg_incident_age <- c(avg_incident_age, mean(incident_age$incident_age))

    # sampling 1:1 case-controls
    interval_surv_data_cases <- interval_data %>%
      filter(incident_disease == 1)
    interval_surv_data_controls <- interval_data %>%
      filter(incident_disease == 0) %>%
      sample_n(size = dim(interval_surv_data_cases)[1] * 1, replace = FALSE)
    incident_1to1_data <- interval_surv_data_cases %>%
      bind_rows(interval_surv_data_controls)
    disease_resid <- resid(lm(data = incident_1to1_data, incident_disease ~ genetic_sex + PC1 + PC2 + PC3 + PC4 + PC5, na.action="na.exclude"))

    for(predict_name_idx in 1:length(predictor_names)){
      predictor <- predictor_names[predict_name_idx]
      incident_joint_1to1_R2 <- cor(pull(incident_1to1_data, !!predictor), disease_resid, use = "complete")^2
      Age_1to2_incident_R2[age_bin_idx,predict_name_idx] <- incident_joint_1to1_R2
      # R2_se_delta_approx <- 2 * sqrt(incident_joint_1to1_R2) * (1-incident_joint_1to1_R2) / sqrt(dim(incident_1to1_data)[1] - 1)
      # se_Age_1to2_incident_R2[age_bin_idx,predict_name_idx] <- R2_se_delta_approx
    }
  }
  colnames(Age_1to2_incident_R2) <- predictor_names
  Age_1to2_incident_R2 <- data.table(age_bin_idx = 1:5, age_median_quintile, avg_incident_age, Age_1to2_incident_R2)
  return(list(prediction_R2, Age_1to2_incident_R2))
}
example_data_incident_only <- example_data %>%
  mutate(disease_age = age_diag) %>%
  mutate(prevelant_disease = if_else(disease_age < age_risk_measure, 1, 0),
         incident_disease = if_else(disease_age >= age_risk_measure, 1, 0)) %>%
  mutate(prevelant_disease = if_else(is.na(prevelant_disease), 0, prevelant_disease),
         incident_disease = if_else(is.na(incident_disease), 0, incident_disease)) %>%
  rename(age_decimal = age_risk_measure)
results <- compute_incident_1to1_manyPredictor(example_data_incident_only, c("RiskScore", "PC1") )
