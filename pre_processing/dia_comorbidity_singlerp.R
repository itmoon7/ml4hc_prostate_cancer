setwd("/Users/dexinli/Dropbox (MIT)/MGH Prostate Research Group/RPDR/Code/Diagnoses/")
library(comorbidity)
library(dplyr)

# Import the data and look at the first six rows
diag_df <- read.csv(file = "diamerged_singlerp.csv")

# Make a new ID variable that represents each patient-day
diag_df$EMPI_day = paste(diag_df$EMPI, diag_df$Date)

# Make subsets of data for each ICD Code_Type
icd9_df <- diag_df[diag_df$ICD9 != -1, ] 
icd10_df <- diag_df[diag_df$Code_Type == "ICD10", ] 

# Calculate comorbidities using charlson index for each subset of data, for each type of code
# this is for each patient-day specifically (bc the package aggregates the diagnoses for each ID)
# and we only want to aggregate within each patient and each day for now
charlson9 <- comorbidity(x = icd9_df, id = "EMPI_day", code = "ICD9", score = "charlson", icd = "icd9", assign0 = FALSE)
charlson10 <- comorbidity(x = icd10_df, id = "EMPI_day", code = "Code", score = "charlson", icd = "icd10", assign0 = FALSE)

# Append the charlson 9 and charlson10 datasets together
charlson_df <- rbind(charlson9, charlson10)

# export it, so I can reimport in pandas and clean there
write.csv(charlson_df,"charlson_singlerp.csv")



