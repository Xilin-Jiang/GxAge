## Age-dependent diseases architecture

This repository provides code and summary data for analysing age-dependent disease architecture, focusing on the genetic (**G**) and environmental (**E**) contributions to complex traits and diseases. It includes simulations, quantitative trait analyses, and binary trait analyses, along with example functions to estimate age-dependent prediction accuracy.

Normally, I would write an R package. However, this project contains multiple parts that are not focused on methodology, so we decided it is better to provide all pipelines as scripts for readability. We do provide a simple [function wrapper](#Estimating-age-dependent-prediction-accuracy) with example data for estimating the age variation of prediction accuracy for a risk score. See below for details.

### Scripts overview

| Script | Purpose | Data requirement |
|--------|---------|------------------|
| `simulations_GxAge.R` | Replicates all simulations | None |
| `QuantitativeTraitsAnalaysis.R` | Estimates age-dependent genetic & environmental variance from summary-level data | Summary-level data only |
| `Binary_traits.R` | Estimates age-dependent prediction accuracy for binary traits | Summary-level data only |

---

## What is age-dependent disease architecture?
Our analysis focuses on both the genetic (**G**) and environmental (**E**) components of complex traits and diseases. 
Previous analyses have focused on **GxAge**, while here we analyse both **GxAge** and **ExAge**. In fact, the main
mechanism we find is ***Exposure Accumulation***, which would not even be detected as **GxAge** in previous analyses (e.g. Robinson et al. 2017 *Nature Genetics*).

We enjoyed reading Miao et al. 2025 *Nature Human Behavior*, which provides a comprehensive description of different **GxE** models (Table 1 of Miao et al.); most of these 
models can be extended by setting **Age = E**. Therefore, we focus on mechanisms that are unique to **Age** and list different models of **GxAge** only as secondary analyses. 

To illustrate why analysing **G** and **E** components of complex traits is important, we use the following schematic figure of the liability-threshold model:

![Schematic Figure of age-dependent disease architecture](Schematic.png)

*Figure: Liability-threshold model showing how genetic (**G**) and environmental (**E**) variance contribute to disease risk as age increases.*

In this figure, disease liability is the sum of the **G** and **E** components. As age increases, more individuals develop the disease by crossing the threshold.
- Both **G** and **E** variances can change with age. Even if only **E** variance increases with age while **G** variance remains constant, the heritability of disease liability will change with age.
- For disease prediction, removing *prevalent cases* will create a negative correlation between **G** and **E**, which reduces the prediction power of both **G** and **E**. 

In this work, we first estimate age-dependent profiles of **G** and **E** variance to show that environmental exposure generally accumulates with age. Then we show
how this impacts the prediction accuracy of both genetic and non-genetic predictors. 

---

## Estimating age-dependent prediction accuracy
After cloning this repository, or simply downloading the `GxAge_functions.R` file if you do not need the example data, you can run the following code to estimate age-dependent prediction accuracy.

The function does five things:

1. Divide the cases into five age quintiles.
2. Create a case-control set for each age quintile so that the case-control ratio and number of 
   samples are the same for each quintile.
3. Regress out age, sex, and covariates from the RiskScore. Then estimate the squared correlation for each quintile.
4. Fit a model of linear age effect on the squared correlations (using the median age of each quintile) and test against a model with constant effect across quintiles. The test is a likelihood ratio test for nested models.
5. Bootstrap across individuals and repeat the procedure to obtain the standard error of the age-effect estimates.

Try this example first:

```r
source("GxAge_functions.R")
example_data <- fread("Summary_data/example_RiskScore_age_data.txt")
# Using flag_remove_prevalent_cases=T to perform incident case prediction
results <- Prediction_R2_Slope_per10years(example_data, flag_normal_transform_RiskScore = T, 
                                          flag_remove_prevalent_cases = T)
print(paste0("Incident case prediction per 10 year change: ", results[[2]]$incident_R2_1to1_slope_10year,
             " z-score = ", results[[2]]$slope_zscore))

# Using flag_remove_prevalent_cases=F to perform prevalent case association
results <- Prediction_R2_Slope_per10years(example_data, flag_normal_transform_RiskScore = T, 
                                          flag_remove_prevalent_cases = F)
print(paste0("Prevalent case association per 10 year change: ", results[[2]]$incident_R2_1to1_slope_10year,
             " z-score = ", results[[2]]$slope_zscore))
```

### Data format
Your data must follow the same format as `example_data`. Columns include:

- **individual_id**: Unique ID for each individual
- **disease**: 1 if disease is observed, 0 if not
- **age_diag**: Age at disease diagnosis; ignored for rows with disease == 0
- **censor_age**: Last observation of the individual
- **RiskScore**: Risk score of interest
- **age_risk_measure**: Age when the risk score was estimated; will be regressed out from the risk score. If not available, set to 0. 
- **genetic_sex**: Genetic sex of the individual; will be regressed out from the risk score. If not available, set to 0. 
- **PC1:PC5**: Additional covariates that should be regressed out. If fewer covariates are available, set unused variables to 0. 

The three age columns (`age_diag`, `censor_age`, and `age_risk_measure`) represent the age of disease onset, the last observation age of an individual, and the age when the risk score was measured. 
- We recommend using the first age at diagnosis for `age_diag` when multiple diagnoses of the disease are present in the individual. 
- `censor_age` should be the last age information that is available: death age when available, or the age when the dataset was last updated. For example, UK Biobank periodically updates their HES dataset, so users can compute the age at the latest update for individuals who are still alive.
- For PRS, `age_risk_measure` can be 0 (at birth); otherwise, it is usually a good choice to use the baseline visit age when risk factors are measured.
- We recommend removing individuals without any information for `age_diag`, as it does not provide useful age information.
- For missing `censor_age`, we recommend using the last age available when an individual interacts with the health system, if it cannot be computed from update dates.

The analysis can still be performed if a small number of age data points are missing. However, if the majority of individuals in your dataset lack information for any of the three age columns, the dataset may not be suitable for longitudinal analysis. 

---

## Simulations
All simulation code is organised in the script `simulations_GxAge.R`. The code below will run the analysis and generate all the main figures in the paper. The purpose of providing this code is to make it flexible for readers to adapt to their own analyses.

```r
source(simulations_GxAge.R)
```

- In the first section of the script (Figure 4B), we show simulations for the EA model and linear liability-threshold model. We estimate prediction accuracy using 
both incident case prediction and prevalent case association.
- In the second section of the script (Figure 4C), we simulate different age dependencies for predictors of different overall prediction accuracy.

---

## Quantitative traits
Code for quantitative traits and disease liability analyses is organised in the script `QuantitativeTraitsAnalaysis.R`. 
We also release all summary estimates of heritability and genetic correlation in the `Summary_data` folder, in addition to the supplementary tables of the paper. 

```r
source(QuantitativeTraitsAnalaysis.R)
```

Running the script will save all the main figures related to quantitative traits or disease liability. 
We recommend users run different sections of the script depending on their own research goals, as there are many estimates in this script. 

---

## Binary/disease traits
Code for disease analyses is organised in the script `Binary_traits.R`. This script focuses on age-dependent prediction accuracy.  
We also release all summary data in the `Summary_data` folder, along with the supplementary tables of the paper. 

```r
source(Binary_traits.R)
```

Running the script will save all the main figures relating to disease traits. 
The beginning parts provide example data that users can follow to estimate age-dependent prediction accuracy. The last section also contains an example showing how to compare the age dependencies of multiple risk scores.

---

## Citation
If you use this code, please cite: [Jiang…Durvasula]().

---

## Contact & contributions
If you have questions, suggestions, or would like to contribute, please contact: **xj262@medschl.cam.ac.uk**.
