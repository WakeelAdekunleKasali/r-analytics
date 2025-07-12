source("01_load_packages.R")
source("02_load_and_clean_data.R")


# Create missing indicator
cancer_mis <- Upd_cancer %>%
  mutate(missing_indicator = ifelse(rowSums(!is.na(dplyr::select(., Bxdiagnosis:'Proliferation Rate'))) == 0, 1, 0))

# Fill diagnosis variable
cancer_dia <- cancer_mis %>%
  group_by(StudyID) %>%
  mutate(diagnosis = Bxdiagnosis) %>%
  fill(diagnosis, .direction = "down") %>%
  ungroup()

# Create baseline-adjusted covariates
cancer_check <- Upd_cancer %>%
  group_by(StudyID) %>%
  mutate(
    real_age = if_else(Timepoints == 0, Age, Age + Time),
    Sex = if_else(Timepoints == 0, `SEX 1=M, 2=F`, NA),
    Ethnicity = if_else(Timepoints == 0, `Ethnicity 1=White, 2=others`, NA),
    lesion_risk_site = if_else(Timepoints == 0, `Lesion Risk Site`, NA),
    Ever_smoker = if_else(Timepoints == 0, `Ever smoker`, NA),
    Ever_alcohol = if_else(Timepoints == 0, `Ever Alcohol`, NA),
    OtherLesPrs = if_else(Timepoints == 0, `Other D Lesion present`, NA)
  ) %>%
  fill(real_age, Sex, Ethnicity, lesion_risk_site, Ever_smoker, Ever_alcohol, OtherLesPrs, .direction = "downup") %>%
  ungroup()

# Create diagnosis_retest
cancer_dia <- cancer_mis %>%
  group_by(StudyID) %>%
  mutate(diagnosis = Bxdiagnosis) %>%
  fill(diagnosis, .direction = "down") %>%
  mutate(diagnosis_retest = ifelse(is.na(diagnosis), NA, ifelse(is.na(Bxdiagnosis), 0, 1))) %>%
  dplyr::select(Timepoints, StudyID, Type, Lesion_area, Bxdiagnosis, diagnosis, diagnosis_retest, Age, missing_indicator) %>%
  ungroup()

# Merge datasets
new_output <- full_join(cancer_check, cancer_dia)

# Reorder columns
prj_data <- new_output %>%
  relocate(diagnosis, .after = Bxdiagnosis) %>%
  relocate(diagnosis_retest, .after = diagnosis) %>%
  relocate(lesion_risk_site, .after = `Lesion Risk Site`) %>%
  relocate(real_age, .after = Age) %>%
  relocate(Sex, .after = `SEX 1=M, 2=F`) %>%
  relocate(Ethnicity, .after = `Ethnicity 1=White, 2=others`) %>%
  relocate(Ever_smoker, .after = `Ever smoker`) %>%
  relocate(Ever_alcohol, .after = `Ever Alcohol`) %>%
  relocate(OtherLesPrs, .after = `Other D Lesion present`) %>%
  relocate(missing_indicator, .after = OtherLesPrs)

spag_Gender <- prj_data[!is.na(cancer$Lesion_area),]
