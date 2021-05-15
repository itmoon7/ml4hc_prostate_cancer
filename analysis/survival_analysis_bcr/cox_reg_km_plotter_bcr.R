library("survival")
library("survminer")
# install.packages("lfe")
# install.packages("hdm")
# install.packages("glmnet")
# install.packages("sandwich")
# install.packages("randomForest")
library(lfe)
library(hdm)
library(glmnet)
library(sandwich)
library(randomForest)

# use non-factorized version to plot semi-parametric survival probabilities 
setwd("~/Downloads/prostate_cancer_project/codes")
prostate_cox_data <- read.csv(file = '../data/df_cox_data_bcr.csv', row.names = 'EMPI')#, header = TRUE, row.names="EMPI")
prostate_cox_data$psa_prior_to_rp <- (prostate_cox_data$psa_prior_to_rp - mean(prostate_cox_data$psa_prior_to_rp))/sqrt(var(prostate_cox_data$psa_prior_to_rp))
prostate_cox_data$Age.at.RP <- (prostate_cox_data$Age.at.RP - mean(prostate_cox_data$Age.at.RP))/sqrt(var(prostate_cox_data$Age.at.RP))

# drop insignificant features 
drop_cols = c('rp_date', 'max_psa', 'min_psa', 'mean_psa', 'ami_agg', 'chf_agg', 'pvd_agg', 'cevd_agg', 'copd_agg', 'rheumd_agg', 'pud_agg', 'mld_agg', 'diabwc_agg', 'hp_agg', 'rend_agg', 'metacanc_agg', 'aids_agg')
prostate_cox_data = prostate_cox_data[,!(names(prostate_cox_data) %in% drop_cols)]

# change column names to be consistent with the survival package
colnames(prostate_cox_data)[grep("bcr_ind", colnames(prostate_cox_data))] <- 'status'
colnames(prostate_cox_data)[grep("time_to_bcr_in_month", colnames(prostate_cox_data))] <- 'time'

# factorize columns
prostate_cox_data$Ethnicity[prostate_cox_data$White == 1] <- 'White'
prostate_cox_data$Ethnicity[prostate_cox_data$Black == 1] <- 'Black'
prostate_cox_data$Ethnicity[prostate_cox_data$Hispanic == 1] <- 'Hispanic'
prostate_cox_data$Ethnicity[is.na(prostate_cox_data$Ethnicity)] <- 'Other/unknown'
prostate_cox_data$Ethnicity <- factor(c(prostate_cox_data$Ethnicity))
prostate_cox_data$Ethnicity <- relevel(prostate_cox_data$Ethnicity, ref = 'Other/unknown')

prostate_cox_data$ptstage[prostate_cox_data$pt2 == 1] <- 'pt2'
prostate_cox_data$ptstage[prostate_cox_data$pt2a == 1] <- 'pt2'
prostate_cox_data$ptstage[prostate_cox_data$pt2b == 1] <- 'pt2'
prostate_cox_data$ptstage[prostate_cox_data$pt2c == 1] <- 'pt2'
prostate_cox_data$ptstage[prostate_cox_data$pt3a == 1] <- 'pt3a'
prostate_cox_data$ptstage[prostate_cox_data$pt3b == 1] <- 'pt3b'
prostate_cox_data <- prostate_cox_data[complete.cases(prostate_cox_data), ]
prostate_cox_data$ptstage <- factor(c(prostate_cox_data$ptstage))
prostate_cox_data$ptstage <- relevel(prostate_cox_data$ptstage, ref = 'pt2')
# prostate_cox_data$ptstage[is.na(prostate_cox_data$ptstage)] <- 'Other stage'

prostate_cox_data$overall_grade_group <- factor(c(prostate_cox_data$overall_grade_group))
prostate_cox_data$margin <- factor(c(prostate_cox_data$margin))
prostate_cox_data$diab_agg <- factor(c(prostate_cox_data$diab_agg))

# rename columns
names(prostate_cox_data)[names(prostate_cox_data) == "overall_grade_group"] <- "Overall grade group"
names(prostate_cox_data)[names(prostate_cox_data) == "margin"] <- "Surgical margin"
names(prostate_cox_data)[names(prostate_cox_data) == "psa_prior_to_rp"] <- "PSA prior"
names(prostate_cox_data)[names(prostate_cox_data) == "diab_agg"] <- "Diabetes"
names(prostate_cox_data)[names(prostate_cox_data) == "ptstage"] <- "pT stage"
names(prostate_cox_data)[names(prostate_cox_data) == "Age.at.RP"] <- "Age"
# names(prostate_cox_data)[names(prostate_cox_data) == "rp_date"] <- "Surgery year"

str(prostate_cox_data)
# prostate_cox_data$overall_grade_group <- relevel(prostate_cox_data$overall_grade_group, ref = 1)

# create a formula
init_flag <- 0
control_vals = c()
ignore_cols = c('time', 'status', 'pT_stage_combined', 'pt2', 'pt2a', 'pt2b', 'pt2c', 'pt3a', 'pt3b', 'mean_psa', 'max_psa', 'min_psa', 'White', 'Black', 'Hispanic')
for (val in colnames(prostate_cox_data)){
  if (!val %in% ignore_cols){
    if (init_flag == 0){
      form_elems <- val
      if (val != 'treat_genetics_matching'){
        control_vals <- c(control_vals, val)}
      init_flag = 1
    } else{
      form_elems <- paste(form_elems, val, sep = ' + ') 
      if (val != 'treat_genetics_matching'){
        control_vals <- c(control_vals, val)}
    }
  }
}


#cup_data_fit
# fit.coxph <- coxph(as.formula(paste("Surv(time, status)", "~", form_elems)), data = prostate_cox_data)#, ties = 'breslow')
# Age + Ethnicity + `Surgical margin` + `prior PSA` + `prior Diabetes` + `Surgery year` + `Overall grade group` + `pT stage`
fit.coxph <- coxph(formula = Surv(time, status) ~ Age + Ethnicity  + `PSA prior` +`Diabetes` + `Overall grade group` + `pT stage` + `Surgical margin`, data = prostate_cox_data)#, ties = 'breslow')
ggforest(fit.coxph, data = prostate_cox_data)

# test proportional hazard assumption
cox.zph(fit.coxph, transform="rank", terms=TRUE, singledf=FALSE, global=TRUE)

df_coefs <- data.frame(fit.coxph$coefficients)
write.csv(df_coefs, '../data/cox_coef_bcr.csv')

# KM curve - Overall grade group based
prostate_cox_data_cm <- read.csv(file = '../data/df_cox_data_bcr.csv', row.names = 'EMPI')#, header = TRUE, row.names="EMPI")
prostate_cox_data_cm$overall_grade_group <- factor(c(prostate_cox_data_cm$overall_grade_group))

prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2 == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2a == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2b == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2c == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3a == 1] <- 'pt3a'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3b == 1] <- 'pt3b'
prostate_cox_data_cm <- prostate_cox_data_cm[complete.cases(prostate_cox_data_cm), ]

colnames(prostate_cox_data_cm)[grep("bcr_ind", colnames(prostate_cox_data_cm))] <- 'status'
colnames(prostate_cox_data_cm)[grep("time_to_bcr_in_month", colnames(prostate_cox_data_cm))] <- 'time'

table(prostate_cox_data_cm$overall_grade_group)

sfit <- survfit(Surv(time, status)~overall_grade_group, data=prostate_cox_data_cm)
# ggsurvplot(sfit)
ggsurv <- ggsurvplot(sfit, pval = FALSE, conf.int = FALSE, legend.title="",
                     xlim = c(0,120),
                     ylim = c(0.6, 1.00),
                     legend.labs=c("1" ,"2",
                                   "3", "4", "5"),
                     risk.table = TRUE,
                     break.time.by = 30,    
                     xlab="Time in months",
                     # surv.median.line = 'hv',
                     ggtheme = theme_classic(),
                     risk.table.height = 0.35,
                     font.title = c(14),#, "bold", "darkblue"),
                     font.subtitle = c(14),# "bold.italic", "purple"),
                     font.caption = c(14), #"plain", "orange"),
                     font.x = c(14),# "bold.italic", "red"),
                     font.y = c(14), #"bold.italic", "darkred"),
                     font.tickslab = c(14))

# chnage legend font size
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size = 14))

ggsurv$table <- ggrisktable(sfit, data = prostate_cox_data_cm,
                            xlim = c(0, 120),
                            break.time.by = 30,   
                            legend.labs=c("1" ,"2",
                                          "3", "4", "5"),
                            color = "strata",
                            font.tickslab = c(14))

ggsurv

surv_diff <- survdiff(Surv(time, status) ~ overall_grade_group, data = prostate_cox_data_cm)
surv_diff

# anova(fit.coxph)

# KM curve - pT stage based
prostate_cox_data_cm <- read.csv(file = '../data/df_cox_data_bcr.csv', row.names = 'EMPI')#, header = TRUE, row.names="EMPI")
drop_pt_stage = c('pt1', 'pt1a', 'pt1b', 'pt1c', 'pt3c', 'pt4', 'pt3')

prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2 == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2a == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2b == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2c == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3a == 1] <- 'pt3a'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3b == 1] <- 'pt3b'
prostate_cox_data_cm <- prostate_cox_data_cm[complete.cases(prostate_cox_data_cm), ]
prostate_cox_data_cm$ptstage <- factor(c(prostate_cox_data_cm$ptstage))


colnames(prostate_cox_data_cm)[grep("bcr_ind", colnames(prostate_cox_data_cm))] <- 'status'
colnames(prostate_cox_data_cm)[grep("time_to_bcr_in_month", colnames(prostate_cox_data_cm))] <- 'time'

prostate_cox_data_cm = prostate_cox_data_cm[!prostate_cox_data_cm$pT_stage_combined %in% drop_pt_stage, ] 
prostate_cox_data_cm$pT_stage_combined = factor(c(prostate_cox_data_cm$pT_stage_combined))

table(prostate_cox_data_cm$pT_stage_combined)
str(prostate_cox_data_cm)

sfit_v2 <- survfit(Surv(time, status)~ptstage, data=prostate_cox_data_cm)
# ggsurvplot(sfit)
ggsurv <- ggsurvplot(sfit_v2, pval = FALSE, conf.int = FALSE, legend.title="",
                     xlim = c(0,120), 
                     ylim = c(0.6, 1.00),
                     legend.labs=c("pt2","pt3a", "pt3b"),
                     risk.table = TRUE,
                     break.time.by = 30,    
                     xlab="Time in months",
                     # surv.median.line = 'hv',
                     ggtheme = theme_classic(),
                     risk.table.height = 0.35,
                     font.title = c(14),#, "bold", "darkblue"),
                     font.subtitle = c(14),# "bold.italic", "purple"),
                     font.caption = c(14), #"plain", "orange"),
                     font.x = c(14),# "bold.italic", "red"),
                     font.y = c(14), #"bold.italic", "darkred"),
                     font.tickslab = c(14))

# chnage legend font size
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size = 14))

ggsurv$table <- ggrisktable(sfit_v2, data = prostate_cox_data_cm,
                            xlim = c(0, 120),
                            break.time.by = 30,   
                            legend.labs=c("pt2" ,"pt3a", "pt3b"),
                            color = "strata",
                            font.tickslab = c(14))

ggsurv

surv_diff <- survdiff(Surv(time, status) ~ pT_stage_combined, data = prostate_cox_data_cm)
surv_diff


# KM curve - margin based
prostate_cox_data_cm <- read.csv(file = '../data/df_cox_data_bcr.csv', row.names = 'EMPI')#, header = TRUE, row.names="EMPI")

prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2 == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2a == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2b == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2c == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3a == 1] <- 'pt3a'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3b == 1] <- 'pt3b'
prostate_cox_data_cm <- prostate_cox_data_cm[complete.cases(prostate_cox_data_cm), ]

colnames(prostate_cox_data_cm)[grep("bcr_ind", colnames(prostate_cox_data_cm))] <- 'status'
colnames(prostate_cox_data_cm)[grep("time_to_bcr_in_month", colnames(prostate_cox_data_cm))] <- 'time'

table(prostate_cox_data_cm$margin)
str(prostate_cox_data_cm)

sfit_v2 <- survfit(Surv(time, status)~margin, data=prostate_cox_data_cm)
# ggsurvplot(sfit)
ggsurv <- ggsurvplot(sfit_v2, pval = FALSE, conf.int = FALSE, legend.title="",
                     xlim = c(0,120), 
                     ylim = c(0.6, 1.00),
                     legend.labs=c("0", "1"),
                     risk.table = TRUE,
                     break.time.by = 30,    
                     xlab="Time in months",
                     # surv.median.line = 'hv',
                     ggtheme = theme_classic(),
                     risk.table.height = 0.25,
                     font.title = c(14),#, "bold", "darkblue"),
                     font.subtitle = c(14),# "bold.italic", "purple"),
                     font.caption = c(14), #"plain", "orange"),
                     font.x = c(14),# "bold.italic", "red"),
                     font.y = c(14), #"bold.italic", "darkred"),
                     font.tickslab = c(14))

# chnage legend font size
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size = 14))

ggsurv$table <- ggrisktable(sfit_v2, data = prostate_cox_data_cm,
                            xlim = c(0, 120),
                            break.time.by = 30,   
                            legend.labs=c("0", "1"),
                            color = "strata",
                            font.tickslab = c(14))

ggsurv

surv_diff <- survdiff(Surv(time, status) ~ margin, data = prostate_cox_data_cm)
surv_diff


# KM curve - diab_agg based
prostate_cox_data_cm <- read.csv(file = '../data/df_cox_data_bcr.csv', row.names = 'EMPI')#, header = TRUE, row.names="EMPI")

prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2 == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2a == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2b == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt2c == 1] <- 'pt2'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3a == 1] <- 'pt3a'
prostate_cox_data_cm$ptstage[prostate_cox_data_cm$pt3b == 1] <- 'pt3b'
prostate_cox_data_cm <- prostate_cox_data_cm[complete.cases(prostate_cox_data_cm), ]

colnames(prostate_cox_data_cm)[grep("bcr_ind", colnames(prostate_cox_data_cm))] <- 'status'
colnames(prostate_cox_data_cm)[grep("time_to_bcr_in_month", colnames(prostate_cox_data_cm))] <- 'time'

table(prostate_cox_data_cm$diab_agg)
str(prostate_cox_data_cm)

sfit_v2 <- survfit(Surv(time, status)~diab_agg, data=prostate_cox_data_cm)
# ggsurvplot(sfit)
ggsurv <- ggsurvplot(sfit_v2, pval = FALSE, conf.int = FALSE, legend.title="",
                     xlim = c(0,120), 
                     ylim = c(0.6, 1.00),
                     legend.labs=c("0", "1"),
                     risk.table = TRUE,
                     break.time.by = 30,    
                     xlab="Time in months",
                     # surv.median.line = 'hv',
                     ggtheme = theme_classic(),
                     risk.table.height = 0.25,
                     font.title = c(14),#, "bold", "darkblue"),
                     font.subtitle = c(14),# "bold.italic", "purple"),
                     font.caption = c(14), #"plain", "orange"),
                     font.x = c(14),# "bold.italic", "red"),
                     font.y = c(14), #"bold.italic", "darkred"),
                     font.tickslab = c(14))

# chnage legend font size
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size = 14))#, color = "black", face = "bold"))

# https://cran.r-project.org/web/packages/survminer/vignettes/Playing_with_fonts_and_texts.html
ggsurv$table <- ggrisktable(sfit_v2, data = prostate_cox_data_cm,
                            xlim = c(0, 120),
                            break.time.by = 30,   
                            legend.labs=c("0", "1"),
                            color = "strata",
                            font.tickslab = c(14))
ggsurv                        
surv_diff <- survdiff(Surv(time, status) ~ diab_agg, data = prostate_cox_data_cm)
surv_diff



# KM curve - unsupervised learning cluster based
prostate_cox_data_cm <- read.csv(file = '../data/df_cox_data_clusters.csv', row.names = 'EMPI')#, header = TRUE, row.names="EMPI")

colnames(prostate_cox_data_cm)[grep("bcr_ind", colnames(prostate_cox_data_cm))] <- 'status'
colnames(prostate_cox_data_cm)[grep("time_to_bcr_in_month", colnames(prostate_cox_data_cm))] <- 'time'

table(prostate_cox_data_cm$cluster_membershp)
str(prostate_cox_data_cm)

sfit_v2 <- survfit(Surv(time, status)~cluster_membershp, data=prostate_cox_data_cm)
# ggsurvplot(sfit)
ggsurv <- ggsurvplot(sfit_v2, pval = FALSE, conf.int = FALSE, legend.title="",
                     xlim = c(0,120), 
                     ylim = c(0.6, 1.00),
                     # legend.labs=c("0", "1"),
                     risk.table = TRUE,
                     break.time.by = 30,    
                     xlab="Time in months",
                     # surv.median.line = 'hv',
                     ggtheme = theme_classic(),
                     risk.table.height = 0.25,
                     font.title = c(14),#, "bold", "darkblue"),
                     font.subtitle = c(14),# "bold.italic", "purple"),
                     font.caption = c(14), #"plain", "orange"),
                     font.x = c(14),# "bold.italic", "red"),
                     font.y = c(14), #"bold.italic", "darkred"),
                     font.tickslab = c(14))

# chnage legend font size
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size = 14))#, color = "black", face = "bold"))

# https://cran.r-project.org/web/packages/survminer/vignettes/Playing_with_fonts_and_texts.html
ggsurv$table <- ggrisktable(sfit_v2, data = prostate_cox_data_cm,
                            xlim = c(0, 120),
                            break.time.by = 30,   
                            # legend.labs=c("0", "1"),
                            color = "strata",
                            font.tickslab = c(14))
ggsurv                        
surv_diff <- survdiff(Surv(time, status) ~ cluster_membershp, data = prostate_cox_data_cm)
surv_diff


# define a helper function for k-fold CV
DML2.for.PLM <- function(z, d, y, dreg, yreg, nfold=2) {
  nobs <- nrow(z) #number of observations
  foldid <- rep.int(1:nfold,times = ceiling(nobs/nfold))[sample.int(nobs)] #define folds indices
  I <- split(1:nobs, foldid)  #split observation indices into folds  
  ytil <- dtil <- rep(NA, nobs)
  cat("fold: ")
  for(b in 1:length(I)){
    dfit <- dreg(z[-I[[b]],], d[-I[[b]]]) #take a fold out
    yfit <- yreg(z[-I[[b]],], y[-I[[b]]]) # take a foldt out
    dhat <- predict(dfit, z[I[[b]],], type="response") #predict the left-out fold 
    yhat <- predict(yfit, z[I[[b]],], type="response") #predict the left-out fold  
    dtil[I[[b]]] <- (d[I[[b]]] - dhat) #record residual for the left-out fold
    ytil[I[[b]]] <- (y[I[[b]]] - yhat) #record residial for the left-out fold
    cat(b," ")
  }
  #rfit <- lm(ytil ~ dtil)    #estimate the main parameter by regressing one residual on the other
  data <- data.frame(cbind(ytil, dtil))
  #   print(head(data))
  rfit <- felm(ytil ~ dtil,data=data) 
  coef.est <- coef(rfit)[2]  #extract coefficient
  #HC <- vcovHC(rfit)
  se    <- summary(rfit,robust=T)$coefficients[2,2] #record robust standard error by County
  cat(sprintf("\ncoef (se) = %g (%g)\n", coef.est , se))  #printing output
  return( list(coef.est =coef.est , se=se, dtil=dtil, ytil=ytil, rfit=rfit) ) #save output and residuals 
}

# control_variables = c("female", "black", "othrace", "factor(dep)", "q1", "q2", "q3", "q4", "q5", "q6", "agelt35", "agegt54", "durable", "lusd", "husd", "muld")
Y <- cup_data[which(colnames(cup_data) %in% c('time', 'status'))]
D <- cup_data[which(colnames(cup_data) == 'treat_genetics_matching')] # 
Z <- cup_data[which(colnames(cup_data) %in% control_vals)]

data <- data.frame(cbind(Y, D, Z))
y <- as.matrix(Y)
d <- as.matrix(D)
z <- as.matrix(Z) #model.matrix( ~.^2, data=Z) #as.matrix(Z)

fit <- glmnet(z, y, family = "cox")

dreg <- function(z,d){ cv.glmnet(z,d,family="gaussian") } #ML method = lasso from glmnet 
yreg <- function(z,y){ cv.glmnet(z,y,family="cox") }  #ML method = lasso from glmnet 
fit.glmnet <- yreg(z, y)
fit.glmnet$beta
DML2.lasso.cv = DML2.for.PLM(z, d, y, dreg, yreg, nfold=10)
