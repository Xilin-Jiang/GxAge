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
