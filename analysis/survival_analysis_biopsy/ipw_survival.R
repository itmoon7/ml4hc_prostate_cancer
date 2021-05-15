# install.packages("IPTWsurvival")
library(ipw)
library(survival)
library("survival")
library("survminer")

setwd("~/Downloads/prostate_cancer_project/codes")
prostate_cox_data <- read.csv(file = '../data/df_cox_data_death_causal_inference_with_sw.csv')#, header = TRUE, row.names="EMPI")

#
prostate_cox_data$ID <- 1:nrow(prostate_cox_data)

# fit KM curve
fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data)
ggsurvplot(fit, data = prostate_cox_data, xlab="Months of follow-up",
           ylab="Survival probability", palette = c("red", "blue"),
           legend.labs=c("Active\nSurveillance", "RP/Radiation"),
           risk.table = TRUE)
# perform log-rank test
surv_diff <- survdiff(Surv(survtime, death) ~ treated, data = prostate_cox_data)
surv_diff

# get survival curves for both treated and untreated
prostate_cox_data_to_plot <- prostate_cox_data[prostate_cox_data$treated == 1, ]
fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)

fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_treated.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_treated.csv", sep = ",")

prostate_cox_data_to_plot <- prostate_cox_data[prostate_cox_data$treated == 0, ]
fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)

fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_untreated.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_untreated.csv", sep = ",")

# factorize columns
prostate_cox_data$Ethnicity[prostate_cox_data$White == 1] <- 'White'
prostate_cox_data$Ethnicity[prostate_cox_data$Black == 1] <- 'Black'
prostate_cox_data$Ethnicity[prostate_cox_data$Asian == 1] <- 'Asian'
prostate_cox_data$Ethnicity[prostate_cox_data$Hispanic == 1] <- 'Hispanic'
prostate_cox_data$Ethnicity[is.na(prostate_cox_data$Ethnicity)] <- 'Other/unknown'
prostate_cox_data$Ethnicity <- factor(c(prostate_cox_data$Ethnicity))
prostate_cox_data$Ethnicity <- relevel(prostate_cox_data$Ethnicity, ref = 'Other/unknown')
table(prostate_cox_data$Ethnicity)
str(prostate_cox_data)

prostate_cox_data$overall_grade_merged <- factor(c(prostate_cox_data$overall_grade_merged))
prostate_cox_data$treated <- factor(c(prostate_cox_data$treated))
prostate_cox_data$ami_agg <- factor(c(prostate_cox_data$ami_agg))
prostate_cox_data$chf_agg <- factor(c(prostate_cox_data$chf_agg))
prostate_cox_data$pvd_agg <- factor(c(prostate_cox_data$pvd_agg))
prostate_cox_data$copd_agg <- factor(c(prostate_cox_data$copd_agg))
prostate_cox_data$rend_agg <- factor(c(prostate_cox_data$rend_agg))
prostate_cox_data$diab_agg <- factor(c(prostate_cox_data$diab_agg))

# rename columns
names(prostate_cox_data)[names(prostate_cox_data) == "overall_grade_merged"] <- "Overall grade group"
names(prostate_cox_data)[names(prostate_cox_data) == "num_pos_cores_sum"] <- "Num. of positive cores"
names(prostate_cox_data)[names(prostate_cox_data) == "auxiiliary_mci_score"] <- "MCI"
# prior comorbidities
names(prostate_cox_data)[names(prostate_cox_data) == "ami_agg"] <- "AMI"
names(prostate_cox_data)[names(prostate_cox_data) == "chf_agg"] <- "CHF"
names(prostate_cox_data)[names(prostate_cox_data) == "pvd_agg"] <- "PVD"
names(prostate_cox_data)[names(prostate_cox_data) == "copd_agg"] <- "COPD"
names(prostate_cox_data)[names(prostate_cox_data) == "diab_agg"] <- "Diabetes"
names(prostate_cox_data)[names(prostate_cox_data) == "rend_agg"] <- "Renal disease"
names(prostate_cox_data)[names(prostate_cox_data) == "Age.at.RP"] <- "Age"
names(prostate_cox_data)[names(prostate_cox_data) == "psa_prior_to_rp"] <- "PSA prior"
names(prostate_cox_data)[names(prostate_cox_data) == "treated"] <- "RP/Radiation"
# names(prostate_cox_data)[names(prostate_cox_data) == "rp_date"] <- "Surgery year"

# estimation of denominator of ip weights
# p.denom <- glm(treated ~ Ethnicity + Age + overall_grade_merged + num_pos_cores_sum
#                + auxiiliary_mci_score + ami_agg + psa_prior_to_rp
#                + chf_agg + pvd_agg + copd_agg
#                + diab_agg + rend_agg, 
#                data=prostate_cox_data, family=binomial())
# prostate_cox_data$pd.treated <- predict(p.denom, prostate_cox_data, type="response")
# 
# summary(prostate_cox_data[which(prostate_cox_data$treated==1),]$pd.treated)
# summary(prostate_cox_data[which(prostate_cox_data$treated==0),]$pd.treated)
# 
# # estimation of numerator of ip weights
# p.num <- glm(treated ~ 1, data=prostate_cox_data, family=binomial())
# prostate_cox_data$pn.treated <- predict(p.num, prostate_cox_data, type="response")
# 
# summary(prostate_cox_data[which(prostate_cox_data$treated==1),]$pn.treated)
# summary(prostate_cox_data[which(prostate_cox_data$treated==0),]$pn.treated)
# 
# # computation of estimated weights
# prostate_cox_data$sw.a <- ifelse(prostate_cox_data$treated==1, prostate_cox_data$pn.treated/prostate_cox_data$pd.treated,
#                                  (1-prostate_cox_data$pn.treated)/(1-prostate_cox_data$pd.treated))
# # prostate_cox_data$sw.a <- ifelse(prostate_cox_data$treated==1, prostate_cox_data$pn.treated/prostate_cox_data$pd.treated,
# #                                  (1-prostate_cox_data$pn.treated)/(1-prostate_cox_data$pd.treated))
# summary(prostate_cox_data$sw.a)
# str(prostate_cox_data)

# + overall_grade_merged  + num_pos_cores_sum + auxiiliary_mci_score + 
#   rp + ami_agg  + chf_agg + pvd_agg + copd_agg + diab_agg + hp_agg + rend_agg +
#   Ethnicity + psa_prior_to_rp
str(prostate_cox_data)
fit.coxph <- coxph(Surv(survtime, death) ~ Age + Ethnicity + `PSA prior` + `RP/Radiation` +`Num. of positive cores` + MCI + `Overall grade group` +
                     `AMI`  + `CHF` + `PVD` + `COPD` + `Diabetes` + `Renal disease`, 
                   data = prostate_cox_data)#, weight = prostate_cox_data$sw)#, ties = 'breslow')
ggforest(fit.coxph, data = prostate_cox_data)

# no weigting results
## treated vs. untreated entire cohort
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_unweighted.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated_unweighted.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_unweighted.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated_unweighted.csv", sep = ",")

# draw survival curve
continuous_cols = c('Age.at.RP', 'auxiiliary_mci_score', 'psa_prior_to_rp')

treat_df = prostate_cox_data[1:2,]
row.names(treat_df) <- c(1,2)

# drop time and status columns 
drops <- c("survtime","death")
treat_df_filtered = treat_df[ , !(names(treat_df) %in% drops)]

cols = c()
for (val in colnames(treat_df_filtered)){
  if (val != "treated"){
    cols = c(cols, val)
  }
}

treat_df_filtered = treat_df_filtered[c("treated", cols)]

for (col in colnames(treat_df_filtered)){
  if (col %in% continuous_cols){
    treat_df_filtered[col] <- mean(prostate_cox_data[[col]])
    #treat_df_filtered[2,col] <- mean(cup_data[cup_data$treat_genetics_matching 
    #                                         == 1,][[col]])
    print(col)
    flush.console()
  } 
  else {
    treat_df_filtered[col] <- mean(prostate_cox_data[[col]])
  }
}

treat_df_filtered$treated <- c(0,1)
treat_df_filtered$overall_grade_merged <- c(1,1)
treat_df_filtered$Age.at.RP <- c(1.5,1.5)
# treat_df_filtered$sw.a = 1.028453
treat_df_filtered
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = treat_df_filtered)
ggsurvplot(fit, censor = FALSE, conf.int = TRUE,# legend.labs=c("Non-matching", "Matching"),
           ggtheme = theme_minimal())
# treat_df_filtered$predicted_cancer_types <- 1
# treat_df_filtered$metastasis_sites <- 1
# treat_df_filtered$tumor_types <- 1
# treat_df_filtered$metastasis_number <- 1

# Survival curves
# fit.coxph is a pre-trained model
# newdata includes two cohorts to compare in the survival curve

## treated vs. untreated entire cohort
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_weighted.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_weighted.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated.csv", sep = ",")

# overall grade group = 1, age > 70, copd (chronic obstructive pulmonary disease), black population 
## treated vs. untreated : overall grade group 1
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1 & prostate_cox_data$overall_grade_merged == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_weighted_og_1.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated_og_1.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
plot(fit)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_treated_pr_weighted_og_1.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_treated_og_1.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0 & prostate_cox_data$overall_grade_merged == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_weighted_og_1.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated_og_1.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
plot(fit)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_untreated_pr_weighted_og_1.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_untreated_og_1.csv", sep = ",")


## age > 70 (age_70_ = 0.5094730266945321)
## treated vs. untreated : age 70
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1 & prostate_cox_data$Age.at.RP >= 0.5094730266945321, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_weighted_age_70_higher.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated_age_70_higher.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_treated_pr_weighted_age_70_higher.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_treated_age_70_higher.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0 & prostate_cox_data$Age.at.RP >= 0.5094730266945321, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_weighted_age_70_higher.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated_age_70_higher.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_untreated_pr_weighted_age_70_higher.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_untreated_age_70_higher.csv", sep = ",")


## treated vs. untreated : black
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1 & prostate_cox_data$Black == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_weighted_black.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated_black.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_treated_pr_weighted_black.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_treated_black.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0 & prostate_cox_data$Black == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_weighted_black.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated_black.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_untreated_pr_weighted_black.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_untreated_black.csv", sep = ",")

## treated vs. untreated : ami
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1 & prostate_cox_data$ami_agg == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_weighted_ami.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated_ami.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_treated_pr_weighted_ami.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_treated_ami.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0 & prostate_cox_data$ami_agg == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_weighted_ami.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated_ami.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_untreated_pr_weighted_ami.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_untreated_ami.csv", sep = ",")


## treated vs. untreated : chf
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1 & prostate_cox_data$chf_agg == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_weighted_chf.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated_chf.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_treated_pr_weighted_chf.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_treated_chf.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0 & prostate_cox_data$chf_agg == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_weighted_chf.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated_chf.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_untreated_pr_weighted_chf.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_untreated_chf.csv", sep = ",")

##
## treated vs. untreated : copd
prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 1 & prostate_cox_data$copd_agg == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/surv_curves_treated_pr_weighted_copd.csv", sep = ",")
write.table(time_steps_treated, file = "../data/time_steps_treated_copd.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_treated_pr_weighted_copd.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_treated_copd.csv", sep = ",")

prostate_cox_data_to_plot = prostate_cox_data[prostate_cox_data$treated == 0 & prostate_cox_data$copd_agg == 1, ]
fit <- survfit(fit.coxph, data = prostate_cox_data, newdata = prostate_cox_data_to_plot)
fit.summary <- summary(fit)

surv_curves_untreated <- fit.summary$surv
time_steps_untreated <- fit.summary$time
write.table(surv_curves_untreated, file = "../data/surv_curves_untreated_pr_weighted_copd.csv", sep = ",")
write.table(time_steps_untreated, file = "../data/time_steps_untreated_copd.csv", sep = ",")

fit <- survfit(Surv(survtime, death) ~ treated, data=prostate_cox_data_to_plot)
fit.summary <- summary(fit)
surv_curves_treated <- fit.summary$surv
time_steps_treated <- fit.summary$time
write.table(surv_curves_treated, file = "../data/km_surv_curves_untreated_pr_weighted_copd.csv", sep = ",")
write.table(time_steps_treated, file = "../data/km_time_steps_untreated_copd.csv", sep = ",")

#

# ggsurvplot(fit, censor = FALSE, conf.int = FALSE, #legend.labs=c("Active surveillance", "Radiation/RP"),
#            ggtheme = theme_minimal())

# base_weight <- ipwpoint(exposure = rp, family = "binomial", link="logit", numerator = ~1,  denominator =~v1+v2+v3....vn, data = base_model, trunc=0.05) #truncation of 5% for few extreme weights if needed

