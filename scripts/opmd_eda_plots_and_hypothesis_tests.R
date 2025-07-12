source("01_load_packages.R")
source("02_load_and_clean_data.R")
source("03_feature_engineering.R")
source("04_risk_analysis.R")


# Prepare overall_data
overall_data <- spag_Gender %>% 
  dplyr::select(StudyID, Timepoints, Lesion_area, Type)

# Medians
medians_overall_data <- overall_data %>%
  dplyr::summarize(mLesion_area = median(Lesion_area, na.rm = TRUE))
median_value <- medians_overall_data$mLesion_area

# Plot for all patients
p1 <- ggplot(overall_data, aes(Timepoints, Lesion_area, group = StudyID)) + 
  geom_point(size = 0.7) +
  geom_line(color = "grey", alpha = 1.2) +
  geom_hline(aes(yintercept = median_value), linetype = "dotted", color = "black") +
  annotate("text", x = 0, y = median_value, label = round(median_value, 3), vjust = -1, color = "blue", angle = 45) +
  labs(x = "Time (months)", y = expression("Lesion area in "~(mm)^2)) +
  ggtitle("A") +
  theme_bw()

# Non-Progressors
nprogressor_data <- spag_Gender %>% filter(Type == "Non-Progressors")
median_value_np <- median(nprogressor_data$Lesion_area, na.rm = TRUE)

p2 <- ggplot(nprogressor_data, aes(Timepoints, Lesion_area, group = StudyID, col = Sex)) + 
  geom_point(size = 0.7) +
  geom_line() +
  geom_hline(aes(yintercept = median_value_np), linetype = "dotted", color = "black") +
  annotate("text", x = 0, y = median_value_np, label = round(median_value_np, 3), vjust = -1, color = "blue", angle = 45) +
  ggtitle("B") +
  theme_bw()

# Progressors
progressor_data <- spag_Gender %>% filter(Type == "Progressors")
median_value_p <- median(progressor_data$Lesion_area, na.rm = TRUE)

p3 <- ggplot(progressor_data, aes(Timepoints, Lesion_area, group = StudyID, col = Sex)) + 
  geom_point(size = 0.7) +
  geom_line() +
  ylim(0, 2000) +
  geom_hline(aes(yintercept = median_value_p), linetype = "dotted", color = "black") +
  annotate("text", x = 0, y = median_value_p, label = round(median_value_p, 3), vjust = -1, color = "blue", angle = 45) +
  ggtitle("C") +
  theme_bw()

# Combine plots
y_limits <- c(0, 3000)
p1 <- p1 + coord_cartesian(ylim = y_limits)
p2 <- p2 + coord_cartesian(ylim = y_limits)
p3 <- p3 + coord_cartesian(ylim = y_limits)
p1 + p2 + p3

# --------------------------
# HYPOTHESIS TESTING SECTION
# --------------------------

# Risk site hypothesis test
obs_diff_prop1 <- cancer_exp %>% 
  specify(formula = progressed ~ risk_site, success = "Yes") %>%
  calculate(stat = "diff in props", order = c("Low", "High"))

set.seed(1234)
rs_null_distribution <- cancer_exp %>%
  specify(formula = progressed ~ risk_site, success = "Yes") %>%
  hypothesize(null="independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in props", order = c("Low", "High"))

p_value1 <- rs_null_distribution %>% 
  get_p_value(obs_stat = obs_diff_prop1, direction = "both")

# Smoking hypothesis test
obs_diff_prop2 <- cancer_exp %>%
  specify(progressed ~ smoking, success = "Yes") %>%
  calculate(stat = "diff in props", order = c("Yes", "No"))

set.seed(1234)
sm_null_distribution <- cancer_exp %>%
  specify(progressed ~ smoking, success = "Yes") %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in props", order = c("Yes", "No"))

p_value2 <- sm_null_distribution %>%
  get_p_value(obs_stat = obs_diff_prop2, direction = "right")

# Alcohol hypothesis test
obs_diff_prop3 <- cancer_exp %>%
  specify(progressed ~ alcohol, success = "Yes") %>%
  calculate(stat = "diff in props", order = c("Yes", "No"))

set.seed(3782)
al_null_distribution <- cancer_exp %>%
  specify(progressed ~ alcohol, success = "Yes") %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in props", order = c("Yes", "No"))

p_value3 <- al_null_distribution %>%
  get_p_value(obs_stat = obs_diff_prop3, direction = "right")

# Other dysplasia hypothesis test
obs_diff_prop4 <- cancer_exp %>%
  specify(progressed ~ oth_ls, success = "Yes") %>%
  calculate(stat = "diff in props", order = c("Yes", "No"))

set.seed(3782)
ot_null_distribution <- cancer_exp %>%
  specify(progressed ~ oth_ls, success = "Yes") %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in props", order = c("Yes", "No"))

p_value4 <- ot_null_distribution %>%
  get_p_value(obs_stat = obs_diff_prop4, direction = "both")

# Gender hypothesis test
obs_diff_prop5 <- cancer_exp %>%
  specify(progressed ~ gender, success = "Yes") %>%
  calculate(stat = "diff in props", order = c("Male", "Female"))

set.seed(1234)
gd_null_distribution <- cancer_exp %>%
  specify(progressed ~ gender, success = "Yes") %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in props", order = c("Male", "Female"))

p_value5 <- gd_null_distribution %>%
  get_p_value(obs_stat = obs_diff_prop5, direction = "both")
