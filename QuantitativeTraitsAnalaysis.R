source("GxAge_functions.R")
# Summary_data contains the estimated heritability and genetic correlation
####################################
# test change of h2g w.r.t. age
####################################
h2g_results <- fread( paste0("Summary_data//h2g_results_agebins.txt"))
age_idx_median <- fread("Summary_data/age_median_baseline_visits.csv")
h2g_results <- h2g_results %>%
  left_join(age_idx_median, by = "age_idx1")

trait_list <- unique(h2g_results$quant_name)
trait_slope_per_year <- c()
slope_se <- c()
h2g_intercept <- c()
likelihood_ratio_trait <- c()
p_value <- c()
for(trait_idx in 1:length(trait_list)){
  trait_name <- trait_list[trait_idx]
  h2g_fit <- h2g_results %>%
    filter(quant_name == !!trait_name)
  # y <- h2g_fit$h2g_mean
  # x <- h2g_fit$age_median
  # se <- h2g_fit$h2g_se
  mle_estimates <- mle_gaussian_with_se(y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se)
  trait_slope_per_year <- c(trait_slope_per_year, mle_estimates[2])
  # # test if the log likelihood is accurate (compare to optim function)
  # optim(par=c(0.1,0.1), fn = function(x) -likelihood_gaussian_with_se(a = x[1], b=x[2],y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se))
  # likelihood_gaussian_with_se(a = mle_estimates[2], b=mle_estimates[1],y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se)
  mls_no_slope_estimate <- mle_no_slope_gaussian_with_se(y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se)
  h2g_intercept <- c(h2g_intercept, mls_no_slope_estimate)
  likelihood_ratio <- -2*(likelihood_gaussian_with_se(a = 0, b=mls_no_slope_estimate,y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se) -
                            likelihood_gaussian_with_se(a = mle_estimates[2], b=mle_estimates[1],y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se))
  likelihood_ratio_trait <- c(likelihood_ratio_trait, likelihood_ratio)
  p_value <- c(p_value, 1-pchisq(likelihood_ratio,df=1))
  # monte-carlo samples for slope se estiamte
  mc_samples <- 100
  h2_samples <- mvtnorm::rmvnorm(mc_samples, mean = h2g_fit$h2g_mean, sigma = diag(h2g_fit$h2g_se^2))
  bt_slope_estimate <- sapply(1:mc_samples, function(x) mle_gaussian_with_se(y = h2_samples[x,], x = h2g_fit$age_median, se = h2g_fit$h2g_se )[2])
  slope_se <- c(slope_se, sd(bt_slope_estimate))
}
h2g_change_test <- data.table(trait_list,
                              trait_slope_per_year,
                              slope_se,
                              h2g_intercept,
                              likelihood_ratio_trait,
                              p_value) %>%
  mutate(bonforonni = p_value * length(trait_list), percent_per_10year = trait_slope_per_year * 10 / h2g_intercept)
h2g_change_test %>%
  filter(bonforonni < 0.1)

Plot_labels <- c("BMI", "Height", "Body WHR","FEV1/FVC",
                 "Smoking Status", "HbA1c", "LDLdirect",
                 "Triglycerides", "Glucose", "Cholesterol",
                 "RBC Width", "MCH", "Eosinophil Count", "Systolic BP",
                 "Diastolic BP", "Reticulocyte Count")
trait_list <- c("bmi","height", "body_WHR", "FEV1_FVC",
                "cov_SMOKING_STATUS", "biochemistry_HbA1c_corrected","biochemistry_LDLdirect_corrected",
                "biochemistry_Triglycerides_corrected", "biochemistry_Glucose","biochemistry_Cholesterol_corrected",
                "blood_RBC_DISTRIB_WIDTH", "blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
                "blood_EOSINOPHIL_COUNT", "Systolic_bp_corrected",
                "Diastolic_bp_corrected", "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT")
match_plot_labs <- data.table(trait_list, Plot_labels)

plot_h2g_change <- h2g_change_test
plot_h2g_change <- plot_h2g_change %>%
  left_join(match_plot_labs, by = c("trait_list") )
plot_h2g_change <- plot_h2g_change %>%
  mutate(percent_per_10year_se = slope_se * 10 / h2g_intercept, significant = bonforonni < 0.1) %>%
  arrange(desc(percent_per_10year))
plot_h2g_change <- plot_h2g_change %>%
  mutate(labels = factor(Plot_labels, levels = plot_h2g_change$Plot_labels))
plt_h2g_change <- ggplot(plot_h2g_change) +
  geom_col(aes(y=percent_per_10year, x=labels, fill = significant), alpha = 0.5) +
  geom_pointrange(aes(
    y=percent_per_10year, ymin=percent_per_10year - percent_per_10year_se,
    ymax = percent_per_10year + percent_per_10year_se,  x=labels, color = significant)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = red, "FALES" = grey)) +
  scale_color_manual(values = c("TRUE" = red, "FALES" = grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(x = "", y = expression(paste("10-year ", ~Delta, h^2, ))) +
  scale_y_continuous(labels = scales::percent)
ggsave(paste0("Figures/h2g_changes_quant.pdf"), plt_h2g_change, width = 6, height = 6)

####################################
# compute var(G) and var(E)
####################################
# note the s.e. requires monte-carelo sampling of individual data, so we only show the raw summary estimates.
# scatter plot
joint_scatter_plot <- list()
joint_scatter_plot[[1]] <- fread(file = paste0("Summary_data/Male_scatterplot_wide_format_G_E_variance_20250127.txt"), sep = "\t")
joint_scatter_plot[[2]] <- fread(file = paste0("Summary_data/Female_scatterplot_wide_format_G_E_variance_20250127.txt"), sep = "\t")
joint_scatter_plot <- bind_rows(joint_scatter_plot)
joint_scatter_plot <- joint_scatter_plot %>%
  left_join(select(h2g_change_test, trait_list, bonforonni)) %>%
  mutate(significant = if_else(bonforonni < 0.05, "H2g_sginificant", "not_significant")) %>%
  mutate(significant = if_else((bonforonni >= 0.1) & (G_bonforonni < 0.1), "G_significant", significant)) %>%
  left_join(match_plot_labs, by = "trait_list") %>%
  mutate(point_label = if_else(bonforonni < 0.1, paste0(Plot_labels," (", sex, ")" ), ""))
plt_G_E_comparison_with_labels_for_reference <- ggplot(joint_scatter_plot) +
  geom_pointrange(aes(
    y=E_percent_per_10year, ymin=E_percent_per_10year - E_percent_per_10year_se,
    ymax = E_percent_per_10year + E_percent_per_10year_se,  x=G_percent_per_10year,
    color = significant, shape = sex),
    size = 0.5) +
  geom_pointrange(aes(
    y=E_percent_per_10year, x=G_percent_per_10year,
    xmin = G_percent_per_10year - G_percent_per_10year_se, xmax = G_percent_per_10year + G_percent_per_10year_se,
    color = significant,  shape = sex),
    size = 0.5) +
  # scale_fill_manual(values = c("H2g_sginificant" = red, "not_significant" = grey, "G_significant" = grey)) +
  scale_color_manual(values = c("H2g_sginificant" = red, "not_significant" = grey, "G_significant" = grey)) +
  scale_shape_manual(values = c("Male" = 0, "Female" = 2)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 25),
        axis.text = element_text(size = 25),
        axis.title= element_text(size = 25) ) +
  labs(x = expression(paste(~Delta, var(G), " w.r.t. age (per 10 year)")), y = expression(paste(~Delta, var(E), " w.r.t. age (per 10 year)"))) +
  geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent) +
  ggrepel::geom_label_repel(data = filter(joint_scatter_plot, point_label != ""),
                            aes(y=E_percent_per_10year, x=G_percent_per_10year, label = point_label),
                            size = 5
  )
ggsave(paste0("Figures/Scatter_with_labels_G_E_variance_diff_changes_quant.pdf"), plt_G_E_comparison_with_labels_for_reference, width = 10, height = 10)

##############################################################
# genetic correlation and phenotypic correlation
##############################################################
rg_results <- fread( paste0("Summary_data/rg_across_agebins.txt")) %>%
  left_join(age_idx_median, by= "age_idx1") %>%
  rename(age_median_1 = age_median ) %>%
  left_join(age_idx_median, by= c("age_idx2" ="age_idx1")) %>%
  rename(age_median_2 = age_median ) %>%
  mutate(age_diff = age_median_2 - age_median_1)

rg_per10year_analyse <- rg_results %>%
  filter(age_idx1 + 1 == age_idx2) %>% # only using 4 pairs
  mutate(rg_per_10_year = exp(log(rg_mean)/(age_median_2 - age_median_1) * 10)) %>%
  mutate(rg_per_10_year_se = rg_per_10_year/rg_mean * rg_se) %>%
  mutate(ivw_weight = 1/rg_per_10_year_se^2) %>%
  group_by(quant_name) %>%
  summarise(rg_per_10_year_mean = sum(rg_per_10_year * ivw_weight)/sum(ivw_weight),
            rg_per_10_year_se = sqrt(1/sum(ivw_weight))) %>%
  arrange(rg_per_10_year_mean) %>%
  left_join(match_plot_labs, by = c("quant_name"="trait_list") )

plot_rg_results <- rg_per10year_analyse %>%
  mutate(rg_diff_p = pnorm((rg_per_10_year_mean-1)/rg_per_10_year_se)) %>%
  mutate(Plot_labels = factor(Plot_labels, levels = rg_per10year_analyse$Plot_labels),
         significant = rg_diff_p * length(unique(rg_results$quant_name)) < 0.1)

plt_rg <- ggplot(plot_rg_results) +
  geom_col(aes(y=rg_per_10_year_mean, x=Plot_labels, fill = significant), alpha = 0.5) +
  geom_pointrange(aes(
    y=rg_per_10_year_mean, ymin=rg_per_10_year_mean - rg_per_10_year_se,
    ymax = rg_per_10_year_mean + rg_per_10_year_se,  x=Plot_labels, color = significant)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = red, "FALES" = grey)) +
  scale_color_manual(values = c("TRUE" = red, "FALES" = grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  # labs(x = "", y = expression(paste("10-year E[",r[g], "]"))) +
  labs(x = "", y = expression(paste("10-year ",rho[g], ))) +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(color = red, size = 1, yintercept = 1, linetype = "dashed")
ggsave(paste0("Figures/per10year_rg_quant.pdf"), plt_rg, width = 6, height = 6)

quant_visit0_1_pheno_cor <- fread(sep = "\t", "Summary_data/Quant_phenotypic_correlation_20240604.txt")

joint_scatter_plot <- quant_visit0_1_pheno_cor %>%
  select(quant_name, per_10_year_cor, se_per_10_year_cor) %>%
  left_join(plot_rg_results, by = "quant_name") %>%
  mutate(rl_diff_p = pnorm((per_10_year_cor-1)/se_per_10_year_cor)) %>%
  mutate(significant = rl_diff_p * length(unique(rg_results$quant_name)) < 0.1) %>%
  mutate(point_label = if_else(significant, paste0(Plot_labels), ""))

plt_rho_G_E_comparison_with_labels_for_reference <- ggplot(joint_scatter_plot) +
  geom_pointrange(data = filter(joint_scatter_plot, significant == T),aes(
    y=rg_per_10_year_mean, ymin=rg_per_10_year_mean - rg_per_10_year_se,
    ymax = rg_per_10_year_mean + rg_per_10_year_se,  x=per_10_year_cor),
    size = 0.5,
    color = red) +
  geom_pointrange(data = filter(joint_scatter_plot, significant == T), aes(
    y=rg_per_10_year_mean, x=per_10_year_cor,
    xmin = per_10_year_cor - se_per_10_year_cor, xmax = per_10_year_cor + se_per_10_year_cor),
    size = 0.5,
    color = red) +
  geom_pointrange(data = filter(joint_scatter_plot, significant == F), aes(
    y=rg_per_10_year_mean, x=per_10_year_cor,
    xmin = per_10_year_cor, xmax = per_10_year_cor),
    size = 0.5,color = grey) +
  # scale_fill_manual(values = c("H2g_sginificant" = red, "not_significant" = grey, "G_significant" = grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  # labs(x = expression(paste("10-year E[",rho[l], "]")), y = expression(paste("10-year E[",r[g], "]"))) +
  labs(x = expression(paste("10-year ",tilde(rho[l]))), y = expression(paste("10-year ",rho[g]))) +
  geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(0.5,1.2))+
  scale_x_continuous(labels = scales::percent, limits = c(0.5,1.2)) +
  ggrepel::geom_label_repel(data = filter(joint_scatter_plot, point_label != ""),
                            aes(y=rg_per_10_year_mean, x=per_10_year_cor, label = point_label)
  )
ggsave(paste0("Figures/Scatter_correlation_10_year_quant.pdf"), plt_rho_G_E_comparison_with_labels_for_reference, width = 6, height = 6)


#########################################################
# h2g and genetic correlation for predicted liability
#########################################################
pop_h2g_predictedliability_results <- fread( paste0("Summary_data/predLiab_h2g_results_all_20250311.txt"))

# load age-bin specific results
h2g_predictedliability_results <- fread( paste0("Summary_data/predLiab_h2g_results_agebins_20241218.txt"))

age_idx_median <- fread("Summary_data/age_median_baseline_visits.csv")

h2g_predictedliability_results <- h2g_predictedliability_results %>%
  left_join(age_idx_median, by = "age_idx1")

trait_list <- unique(h2g_predictedliability_results$quant_name)
trait_slope_per_year <- c()
slope_se <- c()
h2g_intercept <- c()
likelihood_ratio_trait <- c()
p_value <- c()
for(trait_idx in 1:length(trait_list)){
  trait_name <- trait_list[trait_idx]
  h2g_fit <- h2g_predictedliability_results %>%
    filter(quant_name == !!trait_name)
  # y <- h2g_fit$h2g_mean
  # x <- h2g_fit$age_median
  # se <- h2g_fit$h2g_se
  mle_estimates <- mle_gaussian_with_se(y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se)
  trait_slope_per_year <- c(trait_slope_per_year, mle_estimates[2])
  # # test if the log likelihood is accurate (compare to optim function)
  # optim(par=c(0.1,0.1), fn = function(x) -likelihood_gaussian_with_se(a = x[1], b=x[2],y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se))
  # likelihood_gaussian_with_se(a = mle_estimates[2], b=mle_estimates[1],y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se)
  mls_no_slope_estimate <- mle_no_slope_gaussian_with_se(y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se)
  h2g_intercept <- c(h2g_intercept, mls_no_slope_estimate)
  likelihood_ratio <- -2*(likelihood_gaussian_with_se(a = 0, b=mls_no_slope_estimate,y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se) -
                            likelihood_gaussian_with_se(a = mle_estimates[2], b=mle_estimates[1],y = h2g_fit$h2g_mean, x = h2g_fit$age_median, se = h2g_fit$h2g_se))
  likelihood_ratio_trait <- c(likelihood_ratio_trait, likelihood_ratio)
  p_value <- c(p_value, 1-pchisq(likelihood_ratio,df=1))
  # monte-carlo samples for slope se estiamte
  mc_samples <- 100
  h2_samples <- mvtnorm::rmvnorm(mc_samples, mean = h2g_fit$h2g_mean, sigma = diag(h2g_fit$h2g_se^2))
  bt_slope_estimate <- sapply(1:mc_samples, function(x) mle_gaussian_with_se(y = h2_samples[x,], x = h2g_fit$age_median, se = h2g_fit$h2g_se )[2])
  slope_se <- c(slope_se, sd(bt_slope_estimate))
}
h2g_predictliability_change_test <- data.table(trait_list,
                                               trait_slope_per_year,
                                               slope_se,
                                               h2g_intercept,
                                               likelihood_ratio_trait,
                                               p_value) %>%
  mutate(bonforonni = p_value * length(trait_list), percent_per_10year = trait_slope_per_year * 10 / h2g_intercept)
h2g_predictliability_change_test %>%
  filter(bonforonni < 0.05)

h2g_predictliability_change_test <- h2g_predictliability_change_test %>%
  mutate(trait_list = gsub("_"," ", trait_list))

plot_h2g_change <- h2g_predictliability_change_test
plot_h2g_change <- plot_h2g_change %>%
  mutate(percent_per_10year_se = slope_se * 10 / h2g_intercept, significant = bonforonni < 0.1) %>%
  arrange(percent_per_10year)
plot_h2g_change <- plot_h2g_change %>%
  mutate(labels = factor(trait_list, levels = plot_h2g_change$trait_list))
plot_liability_h2g <- ggplot(plot_h2g_change) +
  geom_col(aes(y=percent_per_10year, x=labels, fill = significant), alpha = 0.5) +
  geom_pointrange(aes(
    y=percent_per_10year, ymin=percent_per_10year - percent_per_10year_se,
    ymax = percent_per_10year + percent_per_10year_se,  x=labels, color = significant)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = red, "FALSE" = grey)) +
  scale_color_manual(values = c("TRUE" = red, "FALSE" = grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(x = "", y = expression(paste(~Delta, h^2, " w.r.t. age (per 10 years)"))) +
  scale_y_continuous(labels = scales::percent)

ggsave(paste0("Figures/PredictedLiability_h2g_changes_quant.pdf"), plot_liability_h2g, width = 5, height = 6)


# secondary analysis: using unconstrained LDSC
noConstraintsIntercept_rg_results <- fread(sep = "\t", paste0("Summary_data/nonInterceptConstraints_predLiab_rg_across_agebins_20250312.txt")) %>%
  left_join(age_idx_median, by= "age_idx1") %>%
  rename(age_median_1 = age_median ) %>%
  left_join(age_idx_median, by= c("age_idx2" ="age_idx1")) %>%
  rename(age_median_2 = age_median ) %>%
  mutate(age_diff = age_median_2 - age_median_1)
# normalise rg
rg_per10year_analyse <- noConstraintsIntercept_rg_results %>%
  filter(age_idx1 + 1 == age_idx2) %>%
  mutate(rg_per_10_year = exp(log(rg_mean)/(age_median_2 - age_median_1) * 10)) %>%
  mutate(rg_per_10_year_se = rg_per_10_year/rg_mean * rg_se) %>%
  mutate(ivw_weight = 1/rg_per_10_year_se^2) %>%
  group_by(quant_name) %>%
  summarise(ivw_mean = sum(rg_per_10_year * ivw_weight)/sum(ivw_weight),
            ivw_se = sqrt(1/sum(ivw_weight))) %>%
  rename(rg_per_10_year_mean = ivw_mean,
         rg_per_10_year_se = ivw_se) %>%
  arrange(rg_per_10_year_mean)  %>%
  mutate(quant_name = gsub("_", " ", quant_name))

plot_rg_results <- rg_per10year_analyse %>%
  mutate(rg_diff_p = 2*(1-pnorm( abs(rg_per_10_year_mean-1)/rg_per_10_year_se ) )) %>%
  mutate(quant_name = factor(quant_name, levels = rg_per10year_analyse$quant_name),
         significant = rg_diff_p * length(unique(noConstraintsIntercept_rg_results$quant_name)) < 0.05)

plt_rg <- ggplot(plot_rg_results) +
  geom_col(aes(y=rg_per_10_year_mean, x=quant_name, fill = significant), alpha = 0.5) +
  geom_pointrange(aes(
    y=rg_per_10_year_mean, ymin=rg_per_10_year_mean - rg_per_10_year_se,
    ymax = rg_per_10_year_mean + rg_per_10_year_se,  x=quant_name, color = significant)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = red, "FALES" = grey)) +
  scale_color_manual(values = c("TRUE" = red, "FALES" = grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(x = "", y = expression(paste("10-year ",rho[g], ""))) +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(color = red, size = 1, yintercept = 1, linetype = "dashed")
ggsave(paste0("Figures/noConstraintsIntercept_predicteLiability_per10year_rg.pdf"), plt_rg, width = 5, height = 6)



