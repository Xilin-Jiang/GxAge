source("GxAge_functions.R")
dir.create("Figures")
#######################################################################################
# Figure 4B
#######################################################################################

# part 1: EA incident case prediction model

set.seed(19940110)
incident_AUC_age <- matrix(NA, nrow = 10, ncol = 5)
incident_obs_R2_1_1_matching <- matrix(NA, nrow = 10, ncol = 5)
incident_liablity_R2 <- matrix(NA, nrow = 10, ncol = 5)
incident_logOR_age <- matrix(NA, nrow = 10, ncol = 5)
h2g <- 0.3
E_variance <- 0.5
noise_variance <- 1 - h2g -  E_variance
N <- 100000
prevelance <- 0.02
cumulate_E <- 0.05 # total variance added 20%
noise_ratio <- 1
for(rep in 1:10){
  G <- sqrt(h2g) * rnorm(N)
  E <- sqrt(E_variance) * rnorm(N)
  noise <- sqrt(noise_variance) * rnorm(N)
  L <- G + E + noise
  var(L)

  PRS <- G + noise_ratio * rnorm(N)

  prevelant_cases <- c()
  L_risk_set <- L[setdiff(1:N,prevelant_cases)]
  for(idx_prev in 1:5){
    # get the thresold in the risk set
    thre_L <- quantile(L_risk_set, 1-prevelance/(1-prevelance *(idx_prev- 1) ) )
    # getting all new cases
    incident_cases <- setdiff(which(L > thre_L), prevelant_cases)
    risk_set <- setdiff(1:N, prevelant_cases)
    # print(c(length(risk_set), length(incident_cases)))
    # assign Y
    Y_incident <- rep(0, N)
    Y_incident[incident_cases] <- 1
    incident_AUC_age[rep,idx_prev] <- pROC::auc(response = Y_incident[risk_set], predictor = PRS[risk_set])

    case_idx <- incident_cases
    control_idx <- sample(risk_set, length(incident_cases), replace = F)
    incident_obs_R2_1_1_matching[rep,idx_prev] <- cor(Y_incident[c(case_idx, control_idx)], PRS[c(case_idx, control_idx)])^2
    # for the population prevalence of liability-scale R2, use incidence rate (as the prevelant cases had been removed)
    incident_liablity_R2[rep,idx_prev] <- LiabR2_ascertain(K = mean(Y_incident[risk_set]), P=mean(Y_incident[risk_set]), cor(Y_incident[risk_set], PRS[risk_set])^2)
    incident_logOR_age[rep,idx_prev] <- summary(glm(Y_incident[risk_set] ~ PRS[risk_set], family  = "binomial"))$coefficients[2,1]

    # combine all individuals in the prevalent set
    prevelant_cases <- union(prevelant_cases, which(L > thre_L))
    # adding EA
    L <- L + sqrt(cumulate_E) * rnorm(N)
    # update the risk set
    risk_set <- setdiff(1:N, prevelant_cases)
    L_risk_set <- L[risk_set]
  }
}
# create mean estimates
EA_AUC_age_mean <- colMeans(incident_AUC_age)
EA_obs_R2_1_1_matching_mean <- colMeans(incident_obs_R2_1_1_matching)
EA_liablity_R2_mean <- colMeans(incident_liablity_R2)
EA_logOR_age_mean <- colMeans(incident_logOR_age)

# create se estimates
EA_AUC_age_se <-  apply(incident_AUC_age, 2, function(x) sd(x)/sqrt(10))
EA_obs_R2_1_1_matching_se <-  apply(incident_obs_R2_1_1_matching, 2, function(x) sd(x)/sqrt(10))
EA_liablity_R2_se <-  apply(incident_liablity_R2, 2, function(x) sd(x)/sqrt(10))
EA_logOR_age_se <-  apply(incident_logOR_age, 2, function(x) sd(x)/sqrt(10))

# part 2&3: LinearLT model for incident case prediction & prevalent case association

library(pROC)
AUC_age <- matrix(NA, nrow = 10, ncol = 5)
obs_R2_1_1_matching <- matrix(NA, nrow = 10, ncol = 5)
liablity_R2 <- matrix(NA, nrow = 10, ncol = 5)
logOR_age <- matrix(NA, nrow = 10, ncol = 5)
incident_AUC_age <- matrix(NA, nrow = 10, ncol = 5)
incident_obs_R2_1_1_matching <- matrix(NA, nrow = 10, ncol = 5)
incident_liablity_R2 <- matrix(NA, nrow = 10, ncol = 5)
incident_logOR_age <- matrix(NA, nrow = 10, ncol = 5)
h2g <- 0.3
E_variance <- 0.5
noise_variance <- 1 - h2g -  E_variance
N <- 100000
prevelance_thre <- c(0.04, 0.08, 0.12, 0.16, 0.2) * 0.5
noise_ratio <- 1
for(rep in 1:10){
  G <- sqrt(h2g) * rnorm(N)
  E <- sqrt(E_variance) * rnorm(N)
  noise <- sqrt(noise_variance) * rnorm(N)
  L <- G + E + noise
  var(L)

  PRS <- G + noise_ratio * rnorm(N)
  # clinical_risk_factor <- E + noise_ratio * rnorm(N)

  for(idx_prev in 1:length(prevelance_thre)){
    prevelance <- prevelance_thre[idx_prev]
    thre_L <- quantile(L, 1-prevelance)
    Y <- 1* (L > thre_L)

    # prevalent cases association
    AUC_age[rep,idx_prev] <- pROC::auc(response = Y, predictor = PRS)
    case_idx <- which(Y == 1)
    control_idx <- sample(which(Y == 0), length(case_idx), replace = F)
    obs_R2_1_1_matching[rep,idx_prev] <- cor(Y[c(case_idx, control_idx)], PRS[c(case_idx, control_idx)])^2
    liablity_R2[rep,idx_prev] <- LiabR2(prevelance, cor(Y, PRS)^2)
    logOR_age[rep,idx_prev] <- summary(glm(Y ~ PRS, family  = "binomial"))$coefficients[2,1]

    # incident cases prediction
    case_prev <- 0.02 # matching the first group
    thre_baseline <- quantile(L, 1-prevelance + case_prev)
    thre_incident <- quantile(L, 1-prevelance )
    risk_set <- which(L <= thre_baseline)
    Y_incident <- 1* (L > thre_incident) * (L  <= thre_baseline)

    incident_AUC_age[rep,idx_prev] <- pROC::auc(response = Y_incident[risk_set], predictor = PRS[risk_set])
    case_idx <- which(Y_incident == 1)
    control_idx <- sample(which((Y == 0)&(Y_incident == 0)), length(case_idx), replace = F)
    incident_obs_R2_1_1_matching[rep,idx_prev] <- cor(Y_incident[c(case_idx, control_idx)], PRS[c(case_idx, control_idx)])^2
    # for the population prevalence of liability-scale R2, use incidence rate (as the prevelant cases had been removed)
    incident_liablity_R2[rep,idx_prev] <- LiabR2_ascertain(K = case_prev, P=mean(Y_incident[risk_set]), cor(Y_incident[risk_set], PRS[risk_set])^2)
    incident_logOR_age[rep,idx_prev] <- summary(glm(Y_incident[risk_set] ~ PRS[risk_set], family  = "binomial"))$coefficients[2,1]
  }
}

# create mean estimates
AUC_age_mean <- colMeans(AUC_age)
obs_R2_1_1_matching_mean <- colMeans(obs_R2_1_1_matching)
liablity_R2_mean <- colMeans(liablity_R2)
logOR_age_mean <- colMeans(logOR_age)
incident_AUC_age_mean <- colMeans(incident_AUC_age)
incident_obs_R2_1_1_matching_mean <- colMeans(incident_obs_R2_1_1_matching)
incident_liablity_R2_mean <- colMeans(incident_liablity_R2)
incident_logOR_age_mean <- colMeans(incident_logOR_age)

# create se estimates
AUC_age_mean_se <-  apply(AUC_age, 2, function(x) sd(x)/sqrt(10))
obs_R2_1_1_matching_se <-  apply(obs_R2_1_1_matching, 2, function(x) sd(x)/sqrt(10))
liablity_R2_se <-  apply(liablity_R2, 2, function(x) sd(x)/sqrt(10))
logOR_age_se <-  apply(logOR_age, 2, function(x) sd(x)/sqrt(10))
incident_AUC_age_se <-  apply(incident_AUC_age, 2, function(x) sd(x)/sqrt(10))
incident_obs_R2_1_1_matching_se <-  apply(incident_obs_R2_1_1_matching, 2, function(x) sd(x)/sqrt(10))
incident_liablity_R2_se <-  apply(incident_liablity_R2, 2, function(x) sd(x)/sqrt(10))
incident_logOR_age_se <-  apply(incident_logOR_age, 2, function(x) sd(x)/sqrt(10))

plot_df <- data.table(age_bin = 1:5, prevelance_thre,
                      AUC_age_mean,obs_R2_1_1_matching_mean,liablity_R2_mean,logOR_age_mean,
                      incident_AUC_age_mean, incident_obs_R2_1_1_matching_mean, incident_liablity_R2_mean, incident_logOR_age_mean,
                      AUC_age_mean_se, obs_R2_1_1_matching_se, liablity_R2_se, logOR_age_se,
                      incident_AUC_age_se, incident_obs_R2_1_1_matching_se, incident_liablity_R2_se, incident_logOR_age_se,
                      EA_AUC_age_mean, EA_obs_R2_1_1_matching_mean, EA_liablity_R2_mean, EA_logOR_age_mean,
                      EA_AUC_age_se, EA_obs_R2_1_1_matching_se, EA_liablity_R2_se, EA_logOR_age_se) %>%
  mutate(prevelance_thre = factor(prevelance_thre, levels = prevelance_thre))

plt_liablity_R2 <- ggplot(data=plot_df) +
  geom_pointrange(aes(
    y=liablity_R2_mean, ymin=liablity_R2_mean - liablity_R2_se,
    ymax = liablity_R2_mean + liablity_R2_se,  x=age_bin,
    ), color = purple, alpha = 0.8,
    size = 1) +
  geom_pointrange(aes(
    y=incident_liablity_R2_mean, ymin=incident_liablity_R2_mean - incident_liablity_R2_se,
    ymax = incident_liablity_R2_mean + incident_liablity_R2_se,  x=age_bin,
  ), color = red, alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = liablity_R2_mean),linewidth = 1, color =purple, linetype = 2) +
  geom_line(aes(x = age_bin, y = incident_liablity_R2_mean),linewidth = 1, color =red, linetype = 2) +
  # scale_fill_gradient(limits = c(0,0.7), low = "white", high = red) +
  geom_pointrange(aes(
    y=EA_liablity_R2_mean, ymin=EA_liablity_R2_mean - EA_liablity_R2_se,
    ymax = EA_liablity_R2_mean + EA_liablity_R2_se,  x=age_bin,
  ), color = grey, alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = EA_liablity_R2_mean),linewidth = 1, color =grey, linetype = 2) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(y = expression(paste("liability-scale ",  R^2 )), x = expression(paste( "Age quintiles"))) +
  # geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  # geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(0,NA)) # +
  # scale_x_continuous(labels = scales::percent, limits = c(-0.9,0.3)) # +
  # ggrepel::geom_label_repel(data = filter(plt_predictliabiliy),
  #                           aes(y=PRS_R2_percent_per_10year, x=incident_R2_1to1_slope_10year, label = binary_name)
  # )
ggsave(paste0("Figures/simulation_incidentPredictionVSprevelanceAssociation_liabilityR2_20250320.pdf"), plt_liablity_R2, width = 4, height = 4)

plt_liablity_AUC <- ggplot(data=plot_df) +
  geom_pointrange(aes(
    y=AUC_age_mean, ymin=AUC_age_mean - AUC_age_mean_se,
    ymax = AUC_age_mean + AUC_age_mean_se,  x=age_bin,
  ), color = purple, alpha = 0.8,
  size = 1) +
  geom_pointrange(aes(
    y=incident_AUC_age_mean, ymin=incident_AUC_age_mean - incident_AUC_age_se,
    ymax = incident_AUC_age_mean + incident_AUC_age_se,  x=age_bin,
  ), color = red, alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = AUC_age_mean),linewidth = 1, color =purple, linetype = 2) +
  geom_line(aes(x = age_bin, y = incident_AUC_age_mean),linewidth = 1, color = red, linetype = 2) +
  # scale_fill_gradient(limits = c(0,0.7), low = "white", high = red) +
  geom_pointrange(aes(
    y=EA_AUC_age_mean, ymin=EA_AUC_age_mean - EA_AUC_age_se,
    ymax = EA_AUC_age_mean + EA_AUC_age_se,  x=age_bin,
  ), color = grey, alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = EA_AUC_age_mean),linewidth = 1, color =grey, linetype = 2) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(y = expression(paste("AUC")), x = expression(paste( "Age quintiles"))) +
  # geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  # geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(0.5, NA))

ggsave(paste0("Figures/simulation_incidentPredictionVSprevelanceAssociation_AUC_20250320.pdf"), plt_liablity_AUC, width = 4, height = 4)

plt_R2_1_1 <- ggplot(data=plot_df) +
  geom_pointrange(aes(
    y=obs_R2_1_1_matching_mean, ymin=obs_R2_1_1_matching_mean - obs_R2_1_1_matching_se,
    ymax = obs_R2_1_1_matching_mean + obs_R2_1_1_matching_se,  x=age_bin,
  ), color = purple, alpha = 0.8,
  size = 1) +
  geom_pointrange(aes(
    y=incident_obs_R2_1_1_matching_mean, ymin=incident_obs_R2_1_1_matching_mean - incident_obs_R2_1_1_matching_se,
    ymax = incident_obs_R2_1_1_matching_mean + incident_obs_R2_1_1_matching_se,  x=age_bin,
  ), color = red,alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = obs_R2_1_1_matching_mean),linewidth = 1, color = purple, linetype = 2) +
  geom_line(aes(x = age_bin, y = incident_obs_R2_1_1_matching_mean),linewidth = 1, color = red, linetype = 2) +
  # scale_fill_gradient(limits = c(0,0.7), low = "white", high = red) +
  geom_pointrange(aes(
    y=EA_obs_R2_1_1_matching_mean, ymin=EA_obs_R2_1_1_matching_mean - EA_obs_R2_1_1_matching_se,
    ymax = EA_obs_R2_1_1_matching_mean + EA_obs_R2_1_1_matching_se,  x=age_bin,
  ), color = grey, alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = EA_obs_R2_1_1_matching_mean),linewidth = 1, color =grey, linetype = 2) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(y = expression(paste(R^2, " 1-1 case-control")), x = expression(paste( "Age quintiles"))) +
  # geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  # geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(0,NA))

ggsave(paste0("Figures/simulation_incidentPredictionVSprevelanceAssociation_R21and1matching_20250320.pdf"), plt_R2_1_1, width = 4, height = 4)

plt_logOR <- ggplot(data=plot_df) +
  geom_pointrange(aes(
    y=logOR_age_mean, ymin=logOR_age_mean - logOR_age_se,
    ymax = logOR_age_mean + logOR_age_se,  x=age_bin,
  ), color = purple,alpha = 0.8,
  size = 1) +
  geom_pointrange(aes(
    y=incident_logOR_age_mean, ymin=incident_logOR_age_mean - incident_logOR_age_se,
    ymax = incident_logOR_age_mean + incident_logOR_age_se,  x=age_bin,
  ), color = red ,alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = logOR_age_mean),linewidth = 1, color = purple , linetype = 2) +
  geom_line(aes(x = age_bin, y = incident_logOR_age_mean),linewidth = 1, color = red, linetype = 2) +
  # scale_fill_gradient(limits = c(0,0.7), low = "white", high = red) +
  geom_pointrange(aes(
    y=EA_logOR_age_mean, ymin=EA_logOR_age_mean - EA_logOR_age_se,
    ymax = EA_logOR_age_mean + EA_logOR_age_se,  x=age_bin,
  ), color = grey, alpha = 0.8,
  size = 1) +
  geom_line(aes(x = age_bin, y = EA_logOR_age_mean),linewidth = 1, color =grey, linetype = 2) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  labs(y = expression(paste("log(OR)")), x = expression(paste( "Age quintiles"))) +
  # geom_abline(slope = 1, intercept = 0, color = red, size = 2, alpha = 0.2) +
  # geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(0,NA))

ggsave(paste0("Figures/simulation_incidentPredictionVSprevelanceAssociation_logOR_20250320.pdf"), plt_logOR, width = 4, height = 4)


#######################################################################################
# Figure 4C
#######################################################################################
library(data.table)
h2g_list <- c(0.05, 0.2, 0.5, 0.99)
save_simulation_data <- list()
idx <- 1
rep_number <- 100
for(h2g in h2g_list){
    print(idx)
    #### now make G a predictor, E a predictor (by adding noise)
    R2_1case1control_PRS_ratio <- matrix(NA, nrow = rep_number, ncol = 5)
    logOR_ratio <- matrix(NA, nrow = rep_number, ncol = 5)
    AUC_ratio <- matrix(NA, nrow = rep_number, ncol = 5)
    liablity_R2 <- matrix(NA, nrow = rep_number, ncol = 5)

    # h2g <- 0.3
    E_variance <- 0
    noise_variance <- 1 - h2g -  E_variance
    N <- 100000
    prevelance_thre <- c(0.02, 0.04, 0.06, 0.08, 0.1)
    incident_rate <- 0.01
    # PRS_noise_ratio <- 0.5
    for(rep_idx in 1:rep_number){
      G <- sqrt(h2g) * rnorm(N)
      E <- sqrt(E_variance) * rnorm(N)
      noise <- sqrt(noise_variance) * rnorm(N)
      L <- G + E + noise

      # after collider effect for 3 age groups
      for(age_idx in 1:5){
        prevelance_bin <- prevelance_thre[age_idx]
        # binary 1-1 casse control matching
        case_idx <-  which( (L > quantile(L, 1-incident_rate-prevelance_bin)) &
                              (L < quantile(L, 1-prevelance_bin)) )
        control_idx <-  which(L <= quantile(L, 1-incident_rate-prevelance_bin))
        control_idx <- sample(control_idx, size = length(case_idx), replace = F)

        Y_outcome <- c(rep(1, length(case_idx)), rep(0, length(control_idx))  )
        PRS <- G[c(case_idx, control_idx)] # + ratio_noise * rnorm(N)


        R2_1case1control_PRS_ratio[rep_idx,age_idx] <- cor(PRS, Y_outcome)^2
        logOR_ratio[rep_idx,age_idx] <-  summary(glm(Y_outcome ~ PRS, family  = "binomial"))$coefficients[2,1]
        AUC_ratio[rep_idx,age_idx] <- pROC::auc(response = Y_outcome, predictor = PRS)
        liablity_R2[rep_idx,age_idx] <- LiabR2_ascertain(K = incident_rate, P=mean(Y_outcome), R20 = cor(Y_outcome, PRS)^2)

      }

    }
    # do the same plotting for 1-1 case control R2
    R2_1case1control_PRS_ratio_mean <- colMeans(R2_1case1control_PRS_ratio/R2_1case1control_PRS_ratio[,1])
    R2_1case1control_PRS_ratio_se <- apply(R2_1case1control_PRS_ratio/R2_1case1control_PRS_ratio[,1], 2, function(x) sd(x)/sqrt(rep_number))
    logOR_ratio_mean <- colMeans(logOR_ratio/logOR_ratio[,1])
    logOR_ratio_se <- apply(logOR_ratio/logOR_ratio[,1], 2, function(x) sd(x)/sqrt(rep_number))
    AUC_ratio <- AUC_ratio - 0.5
    AUC_ratio_mean <- colMeans(AUC_ratio/AUC_ratio[,1])
    AUC_ratio_se <- apply(AUC_ratio/AUC_ratio[,1], 2, function(x) sd(x)/sqrt(rep_number))
    liablity_R2_mean <- colMeans(liablity_R2/liablity_R2[,1])
    liablity_R2_se <- apply(liablity_R2/liablity_R2[,1], 2, function(x) sd(x)/sqrt(rep_number))

    save_simulation_data[[idx]] <- data.table(age_bin = 1:5, R2_1case1control_PRS_ratio_mean,
                                       R2_1case1control_PRS_ratio_se,
                                       logOR_ratio_mean,
                                       logOR_ratio_se,
                                       AUC_ratio_mean,
                                       AUC_ratio_se,
                                       liablity_R2_mean,
                                       liablity_R2_se) %>%
      mutate(h2g= h2g)
    print(c(h2g, cor(L[c(case_idx, control_idx)],PRS)^2) )
    idx <- idx + 1
}
save_simulation_data <- bind_rows(save_simulation_data)
df_differentPRS_prediction_data <- save_simulation_data %>%
  mutate(h2g = factor(h2g, c(0.05, 0.2, 0.5, 0.99)))

plt_PRSeffect_differentPRSr2 <- ggplot(df_differentPRS_prediction_data) +
  geom_pointrange(aes(x = age_bin, y = R2_1case1control_PRS_ratio_mean,
                      ymin=R2_1case1control_PRS_ratio_mean - 1.96*R2_1case1control_PRS_ratio_se,
                      ymax=R2_1case1control_PRS_ratio_mean + 1.96*R2_1case1control_PRS_ratio_se,
                      color=h2g), fatten = 0.5, size = 5, alpha = 0.8) +
  geom_line(aes(x = age_bin, y = R2_1case1control_PRS_ratio_mean, color=h2g),linewidth = 1,  linetype = 2) +
  scale_color_manual(values = c(red, blue, green, grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(title = "", x = "Age quintile", y = expression(paste("Relative ", R^2, " (case:control=1:1)")))
ggsave(paste0("Figures/simulation_1to1CaseControl_differentPRSr2.pdf"), plt_PRSeffect_differentPRSr2, width = 4, height = 4)

# AUC
plt_AUC_differentPRSr2 <- ggplot(df_differentPRS_prediction_data) +
  geom_pointrange(aes(x = age_bin, y = AUC_ratio_mean,
                      ymin=AUC_ratio_mean - 1.96*AUC_ratio_se,
                      ymax=AUC_ratio_mean + 1.96*AUC_ratio_se,
                      color=h2g), fatten = 0.5, size = 5, alpha = 0.8) +
  geom_line(aes(x = age_bin, y = AUC_ratio_mean, color=h2g),linewidth = 1,  linetype = 2) +
  scale_color_manual(values = c(red, blue, green, grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(title = "", x = "Age quintile", y = expression(paste("Relative (AUC - 50%)")))
ggsave(paste0("Figures/simulation_1to1AUC_differentPRSr2.pdf"), plt_AUC_differentPRSr2, width = 4, height = 4)

# log OR
plt_logOR_differentPRSr2 <- ggplot(df_differentPRS_prediction_data) +
  geom_pointrange(aes(x = age_bin, y = logOR_ratio_mean,
                      ymin=logOR_ratio_mean - 1.96*logOR_ratio_se,
                      ymax=logOR_ratio_mean + 1.96*logOR_ratio_se,
                      color=h2g), fatten = 0.5, size = 5, alpha = 0.8) +
  geom_line(aes(x = age_bin, y = logOR_ratio_mean, color=h2g),linewidth = 1,  linetype = 2) +
  scale_color_manual(values = c(red, blue, green, grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(title = "", x = "Age quintile", y = expression(paste("Relative log(OR)")))
ggsave(paste0("Figures/simulation_logOR_differentPRSr2.pdf"), plt_logOR_differentPRSr2, width = 4, height = 4)

# liability R2
plt_liabR2_differentPRSr2 <- ggplot(df_differentPRS_prediction_data) +
  geom_pointrange(aes(x = age_bin, y = liablity_R2_mean,
                      ymin=liablity_R2_mean - 1.96*liablity_R2_se,
                      ymax=liablity_R2_mean + 1.96*liablity_R2_se,
                      color=h2g), fatten = 0.5, size = 5, alpha = 0.8) +
  geom_line(aes(x = age_bin, y = liablity_R2_mean, color=h2g),linewidth = 1,  linetype = 2) +
  scale_color_manual(values = c(red, blue, green, grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(title = "", x = "Age quintile", y = expression(paste("Relative liability-scale ", R^2)))
ggsave(paste0("Figures/simulation_liabR2_differentPRSr2.pdf"), plt_liabR2_differentPRSr2, width = 4, height = 4)


plt_PRSeffect_differentPRSr2 <- ggplot(df_differentPRS_prediction_data) +
  geom_pointrange(aes(x = age_bin, y = R2_1case1control_PRS_ratio_mean,
                      ymin=R2_1case1control_PRS_ratio_mean - 1.96*R2_1case1control_PRS_ratio_se,
                      ymax=R2_1case1control_PRS_ratio_mean + 1.96*R2_1case1control_PRS_ratio_se,
                      color=h2g), fatten = 0.5, size = 5, alpha = 0.8) +
  geom_line(aes(x = age_bin, y = R2_1case1control_PRS_ratio_mean, color=h2g),linewidth = 1,  linetype = 2) +
  scale_color_manual(values = c(red, blue, green, grey)) +
  theme(legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15) ) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(title = "", x = "Age quintile", y = expression(paste("Relative ", R^2, " (case:control=1:1)")))

ggsave(paste0("Figures/simulation_1to1CaseControl_differentPRSr2.pdf"), plt_PRSeffect_differentPRSr2, width = 4, height = 4)
