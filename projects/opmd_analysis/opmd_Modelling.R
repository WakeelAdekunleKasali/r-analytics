
# Load required packages
library(dplyr)
library(ggplot2)
library(ggpubr)

# Function to assign timepoint intervals
assign_timepoints <- function(df) {
  df %>%
    mutate(timepoint = case_when(
      Time_point >= 0 & Time_point <= 4 ~ 1,
      Time_point >= 5 & Time_point <= 9 ~ 2,
      Time_point >= 10 & Time_point <= 14 ~ 3,
      Time_point >= 15 & Time_point <= 19 ~ 4,
      Time_point >= 20 & Time_point <= 25 ~ 5,
      Time_point >= 26 & Time_point <= 33 ~ 6,
      Time_point >= 34 & Time_point <= 40 ~ 7,
      Time_point >= 41 & Time_point <= 46 ~ 8,
      Time_point >= 47 & Time_point <= 53 ~ 9,
      Time_point >= 54 & Time_point <= 62 ~ 10,
      Time_point >= 63 & Time_point <= 72 ~ 11,
      Time_point >= 73 & Time_point <= 80 ~ 12,
      Time_point >= 81 & Time_point <= 90 ~ 13,
      Time_point >= 91 & Time_point <= 105 ~ 14,
      Time_point >= 106 & Time_point <= 178 ~ 15,
      TRUE ~ NA_real_
    ))
}

# Function to compute Mann-Whitney U statistic
mann_whitney_u <- function(df) {
  df <- df %>% mutate(Rank = rank(Lesion_area, ties.method = 'average'))
  g1 <- df %>% filter(Type == "Progressors")
  g2 <- df %>% filter(Type == "Non-Progressors")
  R1 <- sum(g1$Rank); R2 <- sum(g2$Rank)
  n1 <- nrow(g1); n2 <- nrow(g2)
  U1 <- n1 * n2 + (n1 * (n1 + 1)) / 2 - R1
  U2 <- n1 * n2 + (n2 * (n2 + 1)) / 2 - R2
  min(U1, U2)
}

# Function to perform permutation test and return p-value
mann_whitney_test_with_permutation <- function(df, n_perm = 5000, seed = 6789) {
  set.seed(seed)
  obs_U <- mann_whitney_u(df)
  perm_Us <- replicate(n_perm, mann_whitney_u(df %>% mutate(Type = sample(Type))))
  pval <- mean(perm_Us <= obs_U)
  return(pval)
}

# Function to generate annotated boxplot
generate_boxplot <- function(df, pval, title) {
  ggplot(df, aes(x = factor(Type), y = Lesion_area)) +
    geom_boxplot() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = title, x = "Group", y = expression("Lesion area in"~ mm^2)) +
    geom_text(
      aes(x = Inf, y = Inf, label = paste("p-value =", format.pval(pval, digits = 3))),
      inherit.aes = FALSE,
      hjust = "inward", vjust = "inward", color = "grey40"
    )
}

# Main analysis block
run_analysis <- function(prj_data_new_al) {
  time_labels <- c(
    "0-4", "5-9", "10-14", "15-19", "20-25", "26-33", "34-40",
    "41-46", "47-53", "54-62", "63-72", "73-80", "81-90",
    "91-105", "106-178"
  )

  study_data <- prj_data_new_al %>%
    select(StudyID, Type) %>%
    distinct()

  prj_data_new_al <- assign_timepoints(prj_data_new_al)

  avg_lesion <- prj_data_new_al %>%
    group_by(StudyID, Type, timepoint) %>%
    summarise(Lesion_area = mean(Lesion_area, na.rm = TRUE), .groups = "drop")

  final_data <- study_data %>%
    left_join(avg_lesion, by = c("StudyID", "Type"))

  plots <- list()
  for (t in 1:15) {
    temp_data <- final_data %>%
      filter(!is.na(Lesion_area), timepoint == t) %>%
      select(Lesion_area, Type)

    if (nrow(temp_data) > 0 && length(unique(temp_data$Type)) == 2) {
      pval <- mann_whitney_test_with_permutation(temp_data)
      plot_title <- paste("Study time", t, "-", time_labels[t], "months")
      plots[[t]] <- generate_boxplot(temp_data, pval, plot_title)
    }
  }

  combined_plot <- ggarrange(plotlist = plots, ncol = 4, nrow = 4)
  return(combined_plot)
}
