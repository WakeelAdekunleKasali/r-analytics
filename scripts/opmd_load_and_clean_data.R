source("01_load_packages.R")


# Read data
orig_data <- read_excel("../Data/New_Dataset.xlsx")
cancer <- read_excel("../Data/New_Dataset.xlsx", sheet = "cancer(cleaned)")
cancer.id <- read_excel("../Data/New_Dataset.xlsx", sheet = "cancer.id")

# Adjust lesion values
cancer$Lesion_area[cancer$StudyID == 2094 & cancer$Timepoints == 25] <- 0
cancer$Lesion_area[cancer$StudyID == 3003 & cancer$Timepoints == 122] <- 0
cancer$Lesion_area[cancer$StudyID == 3011 & cancer$Timepoints == 28] <- 0
Upd_cancer <- cancer

# Remove NAs for plotting
Upd_cancer_spg <- Upd_cancer[!is.na(cancer$Lesion_area),]
