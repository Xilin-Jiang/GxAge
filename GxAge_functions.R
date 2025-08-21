library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(mvtnorm)
library(survival)
library(data.table)
library(dplyr)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]

LiabR2<-function(K, R20){

  #population prevalence
  K=K

  #assume follow Normal distribution
  #the threshold on the normal distribution which truncates the proportion of disease prevalence
  thd=qnorm(1-K)
  zv=dnorm(thd)

  #observed n2
  R2O=R20

  #transfer using equation (8) of Lee
  R2=R2O*K*(1-K)/zv^2
  return(R2)
}

# if the case controls are ascertained, using following model, according to equation 14 of Lee 2012 Genetic epidemiology
LiabR2_ascertain <-function(K, P, R20){
  #K=pop prevalence
  #P=proportion of cases in study
  #R20 observe scale R2
  thd <-  -qnorm(K,0,1)
  zv <- dnorm(qnorm(K))
  mv <- zv/K
  theta <- mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd)
  cv <- (K * ( 1 - K)/zv)^2 / (P * (1-P))
  R2 <- R20 * cv/(1+R20*theta*cv) # note this assumes genetic variance is the same in population vs. case-control, which is generally true
  return(R2)
}

# X <- rnorm(1000)
# Y <- X + rnorm(1000)
# K <- 0.1
# P <- 0.25
# thre_Y <- quantile(Y, 1-K)
# D <- (Y > thre_Y) * 1
# joint_data <- data.table(X, D)
# joint_cases <- filter(joint_data, D == 1)
# joint_controls <- filter(joint_data, D == 0) %>%
#   sample_n(dim(joint_cases)[1]*(1-P)/P)
# joint_data <- joint_cases %>%
#   bind_rows(joint_controls)
# lm_md <- summary(lm(data = joint_data, formula = D ~ X))
# R20 <- lm_md$r.squared
# LiabR2_ascertain(K, P, R20)
# summary(lm(formula = Y ~ X))$r.squared

# computing correlation between E from phenotypic correlation
corr_E_computation <- function(h2g, phe_corr, rg = 1){
  # important simplifying assumption: h2g is the same at the two age points, rg = 1
  cov_G <- h2g * rg
  return((phe_corr - cov_G)/(1-h2g))
}

# likelihood function with slope
likelihood_gaussian_with_se <- function(y, x, a, b, se){ # loglikelihood of y ~ N(ax+b, se^2)
  # sum(log(dnorm(y, mean=a*x + b, sd = se)))
  ll <- sum(-1/(2*se^2) * (y - a*x - b)^2) - 1/2 * sum(log(2 * pi * se^2) )
  return(ll)
}

# fit model y ~ N(ax+b, se^2)
mle_gaussian_with_se <- function(y, x, se){ # loglikelihood of y ~ N(ax+b, se^2)
  sigma_se <-  diag(1/se^2)
  XT <- rbind(rep(1, length(x)), x)
  return(solve(XT %*% sigma_se %*% t(XT)) %*% XT %*% sigma_se %*% y )
}

mle_no_slope_gaussian_with_se <- function(y, x, se){ # loglikelihood of y ~ N(ax+b, se^2)
  return(sum(y/se^2)/sum(1/se^2))
}

# compute the cov(G,E)/h2, which is the proportion of decrease indicated by the prevelance and h2 under liability threshold model
relative_cov_collider <- function(h2, prev){
  if(prev == 0){
    cov_div_h2  <- 0
  }else{
    threshold_gaussian <- qnorm(1-prev)
    cov_div_h2 <- -(1-h2)/sqrt(2*pi) * exp(-threshold_gaussian^2/2) * (threshold_gaussian + exp(-threshold_gaussian^2/2)/sqrt(2*pi))
  }
  return(cov_div_h2)
}

# functions below are for computing slope of prediction accuracy
slope_collider <- function(interval_censored_df){
  case_age <- interval_censored_df %>%
    filter(disease == 1)
  age_quintile <- quantile(case_age$age_decimal, probs = seq(0,1,0.2))
  age_median_quintile <- quantile(case_age$age_decimal, probs = seq(0.1,1,0.2))
  observed_R2 <- c()
  observed_R2_se <- c()
  prevelant_diseae_proportion <- c()
  age_diff_liability_case <- c()
  number_cases_per_bin <- floor(length(case_age$age_decimal)/5)
  for(age_bin_idx in 1:5){
    age_start <- age_quintile[age_bin_idx]
    age_end <- age_quintile[age_bin_idx + 1]
    interval_surv_data <- interval_censored_df %>%
      data.table() %>%
      filter(age_decimal >= age_start) %>% # filter out all who had disease
      dplyr::mutate(disease_status = (age_decimal < age_end)*disease ) %>%
      mutate(observe_age = ifelse(disease_status, age_decimal,  pmin(age_decimal, age_end) ))   %>%
      ungroup()

    # sampling 1:1 case-controls
    interval_surv_data_cases <- interval_surv_data %>%
      filter(disease_status == 1)
    interval_surv_data_controls <- interval_surv_data %>%
      filter(disease_status == 0) %>%
      sample_n(size = number_cases_per_bin * 1, replace = FALSE)
    interval_surv_data <- interval_surv_data_cases %>%
      bind_rows(interval_surv_data_controls)

    # compute the residuals then corr
    disease_resid <- resid(lm(data = interval_surv_data, disease_status ~ genetic_sex  + PC1 + PC2 + PC3 + PC4 + PC5, na.action="na.exclude"))
    cor_predliability_diseaseResid <- cor(disease_resid, interval_surv_data$RiskScore_age_sex_covar_adjusted) # liability_transformed)
    # cor_predliability_diseaseResid <- cor(interval_surv_data$disease_status, interval_surv_data$RiskScore_age_sex_covar_adjusted)
    observed_R2 <- c(observed_R2, cor_predliability_diseaseResid^2)
    R2_se_delta_approx <- 2 * abs(cor_predliability_diseaseResid) * (1-cor_predliability_diseaseResid^2) / sqrt(dim(interval_surv_data)[1] - 1)
    observed_R2_se <- c(observed_R2_se, R2_se_delta_approx)

    # compute cases that are prevelant
    prevelant_diseae_proportion <- c(prevelant_diseae_proportion, mean(interval_surv_data_cases$prevelant_disease))
    age_diff_liability_case <- c(age_diff_liability_case, mean(interval_surv_data_cases$age_decimal - interval_surv_data_cases$age_risk_measure))
  }
  liability_binary_R2_each_age <- data.table(age_bin_idx = 1:5, observed_R2,
                                             observed_R2_se,
                                             prevelant_diseae_proportion,
                                             age_diff_liability_case)

  # observed R2
  R2_slope_estimates <- mle_gaussian_with_se(y = liability_binary_R2_each_age$observed_R2,
                                             x = age_median_quintile,
                                             se = liability_binary_R2_each_age$observed_R2_se)
  R2_slope_per_year <- R2_slope_estimates[2]
  R2_no_slope_estimate <- mle_no_slope_gaussian_with_se(y = liability_binary_R2_each_age$observed_R2,
                                                        x = age_median_quintile,
                                                        se = liability_binary_R2_each_age$observed_R2_se)
  likelihood_ratio <- -2*(likelihood_gaussian_with_se(a = 0, b=R2_no_slope_estimate,
                                                      y = liability_binary_R2_each_age$observed_R2,
                                                      x = age_median_quintile,
                                                      se = liability_binary_R2_each_age$observed_R2_se) -
                            likelihood_gaussian_with_se(a = R2_slope_estimates[2], b=R2_slope_estimates[1],
                                                        y = liability_binary_R2_each_age$observed_R2,
                                                        x = age_median_quintile,
                                                        se = liability_binary_R2_each_age$observed_R2_se))
  R2_p_value <-  1-pchisq(likelihood_ratio,df=1)

  # point estimates
  interval_surv_data <- interval_censored_df %>%
    data.table() %>%
    mutate(disease_status = disease ) %>%
    mutate(observe_age = age_decimal)   %>%
    ungroup()
  # sampling 1:1 case-controls
  interval_surv_data_cases <- interval_surv_data %>%
    filter(disease_status == 1)
  interval_surv_data_controls <- interval_surv_data %>%
    filter(disease_status == 0) %>%
    sample_n(size = dim(interval_surv_data_cases)[1] * 1, replace = FALSE)
  interval_surv_data <- interval_surv_data_cases %>%
    bind_rows(interval_surv_data_controls)


  # compute the residuals then corr
  disease_resid <- resid(lm(data = interval_surv_data, disease_status ~ genetic_sex + PC1 + PC2 + PC3 + PC4 + PC5, na.action="na.exclude"))
  cor_predliability_diseaseResid <- cor(disease_resid, interval_surv_data$RiskScore_age_sex_covar_adjusted) # liability_transformed)
  # cor_predliability_diseaseResid <- cor(interval_surv_data$disease_status, interval_surv_data$RiskScore_age_sex_covar_adjusted)
  all_observed_R2 <- cor_predliability_diseaseResid^2
  all_observed_R2_se <- 2 * abs(cor_predliability_diseaseResid) * (1-cor_predliability_diseaseResid^2) / sqrt(dim(interval_surv_data)[1] - 1)

  slope_results <- data.table(R2_slope_per_year, R2_no_slope_estimate, R2_p_value,
                              all_observed_R2, all_observed_R2_se)
  return(list(liability_binary_R2_each_age, slope_results))
}
Prediction_R2_Slope_per10years <- function(longitidinal_disease_RiskScore_data,
                                           flag_normal_transform_RiskScore = T,
                                           flag_remove_prevalent_cases = T){
  # controlling the covariates
  # Note age of Risk Scores has to be regressed out as under liability-threshold model, we are interested in the
  # non-age dependent component of risk score
  longitidinal_disease_RiskScore_data$RiskScore_age_sex_covar_adjusted <- resid(lm(data = longitidinal_disease_RiskScore_data, RiskScore ~ genetic_sex + age_risk_measure + PC1 + PC2 + PC3 + PC4 + PC5, na.action="na.exclude"))

  if(flag_normal_transform_RiskScore){
    # get a normal distribution of risk scores; not required but this fits the liability-threshold model
    r <- rank(longitidinal_disease_RiskScore_data$RiskScore_age_sex_covar_adjusted, na.last = "keep")
    n <- sum(!is.na(longitidinal_disease_RiskScore_data$RiskScore_age_sex_covar_adjusted))
    # Map the ranks to quantiles of the standard normal distribution
    longitidinal_disease_RiskScore_data$RiskScore_age_sex_covar_adjusted <- qnorm((r - 0.5) / n)
  }


  longitidinal_disease_RiskScore_data <- longitidinal_disease_RiskScore_data %>%
    mutate(disease_age = age_diag) %>%
    mutate(prevelant_disease = if_else(disease_age < age_risk_measure, 1, 0)) %>%
    mutate(prevelant_disease = if_else(is.na(prevelant_disease), 0, prevelant_disease)) %>%
    mutate(disease = replace_na(disease, 0), age_decimal = if_else(is.na(age_diag), censor_age, age_diag)) %>%
    select(- age_diag, -disease_age)

  if(flag_remove_prevalent_cases){
    longitidinal_disease_RiskScore_data <- longitidinal_disease_RiskScore_data %>%
      filter(prevelant_disease == 0 )
  }

  results <- slope_collider(longitidinal_disease_RiskScore_data)
  liability_binary_R2_each_age <- results[[1]]
  slopes_each_trait <- results[[2]]

  mc_samples <- 50
  bt_slopes_each_trait <- list()
  for(bt_idx in 1:mc_samples){
    print(paste0("bootstrapping sample: ", bt_idx))
    bt_interval_censored_df <- longitidinal_disease_RiskScore_data %>%
      slice_sample(prop = 1, replace = T)
    bt_results <- slope_collider(bt_interval_censored_df)
    bt_slopes_each_trait[[bt_idx]] <- bt_results[[2]] %>%
      mutate( bt_idx = bt_idx) %>%
      select(bt_idx, everything())
  }

  results_bt_se_slopes_each_trait <- bind_rows(bt_slopes_each_trait) %>%
    summarise(across(-c("bt_idx"), sd)) %>%
    setNames(paste0('bootstrapping_se_', names(.)))

  Prediction_Slope_estimates <- cbind(slopes_each_trait, select(results_bt_se_slopes_each_trait,
                                                                bootstrapping_se_R2_slope_per_year,
                                                                bootstrapping_se_all_observed_R2)) %>%
    mutate(incident_R2_1to1_slope_10year = R2_slope_per_year * 10 / all_observed_R2,
           incident_R2_1to1_slope_10year_se = bootstrapping_se_R2_slope_per_year * 10 / all_observed_R2) %>%
    mutate(slope_zscore = incident_R2_1to1_slope_10year/incident_R2_1to1_slope_10year_se)
  return(list(liability_binary_R2_each_age, Prediction_Slope_estimates))
}


