source("01_load_packages.R")
source("02_load_and_clean_data.R")
source("03_feature_engineering.R")


# Analysis of categorical variables
cancer.id <- read_excel("../Data/New_Dataset.xlsx", sheet = "cancer.id")
cancer_exp <- cancer.id %>%
  mutate(
    progressed = recode_factor(Progression, `1` = "Yes", `0` = "No"),
    risk_site = recode_factor(`Lesion risk site`, `1` = "High", `0` = "Low"),
    smoking = recode_factor(Smoke, `1` = "Yes", `0` = "No"),
    alcohol = recode_factor(Alcohol, `1` = "Yes", `0` = "No"),
    oth_ls = recode_factor(`Other D Lesion present`, `1` = "Yes", `0` = "No"),
    gender = recode_factor(Gender, `1` = "Yes", `0` = "No")
  )

data1 <- matrix(c(2, 2, 24, 11), 2, 2, byrow = TRUE)
dimnames(data1) <- list("X" = c("Low", "High"), "Y" = c("Non Progressor","Progressor"))
RR <- riskratio(data1, method = "wald")
OR <- oddsratio(data1, method = "wald")
OddsRatioSite <- round(OR$measure[2, 1], 2)

obs_diff_prop1 <- cancer_exp %>%
  specify(progressed ~ risk_site, success = "Yes") %>%
  calculate(stat = "diff in props", order = c("Low", "High"))

set.seed(1234)
rs_null_distribution <- cancer_exp %>%
  specify(progressed ~ risk_site, success = "Yes") %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in props", order = c("Low", "High"))

p_value1 <- rs_null_distribution %>%
  get_p_value(obs_stat = obs_diff_prop1, direction = "both")

percentile_ci_1 <- rs_null_distribution %>%
  get_confidence_interval(level = 0.95, type = "percentile")
