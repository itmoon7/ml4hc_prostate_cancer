# Using Machine Learning to Uncover Prognosis of Patients with Localized Prostate Cancer

Abstract : 
Prostate cancer is one of the most commonly diagnosed cancer among males in the US. Although survival rates of patients diagnosed with localized disease are fairly high following Radical Prostatectomy (RP), about 10 - 15 \% patients experience biochemical recurrence (BCR). To gain clinical insights into the prognosis of patients with the disease, it is crucial to utilize Electronic Health Records (EHR) to build a data based clinical decision support system. To that end, we built the rule-based NLP algorithm which extracts clinical features in the unstructured pathology notes, obtained from Research Patient Data Registry, a centralized clinical data registry that gathers clinical information from various hospital systems. We validated our model's performance with a clinician and the model achieved overall 90\% accuracy on all abstracted features from 20 randomly selected pathology reports. Using those features along with other clinical information including age, ethnicity, and comorbidity, we performed survival analysis on the cohort of prostate cancer patients (n = 1,751) who underwent RP and discovered the extracted clinical features such as overall grade group and pathological tumor stages were significantly associated with the increased risk of BCR. Furthermore, we used the doubly-robust causal inference framework to estimate treatment effect of RP/Radiation on the survival of the cohort of patients with prostate biopsy pathology records demonstrating malignancy (n = 3,945). We estimates heterogeneous group treatment effects across different risk groups defined by age, overall disease stage, and comorbidities and found the strong evidence of survival benefit of the treatment across the groups. 

Team Members :
Dexin Li, Intae Moon, Madhur Nayan, Ashwin Srinivasan

Requires pacakges :
- R : survival, survminer
- Python : pandas, numpy, tqdm, xgboost, sklearn, matplotlib, doubleml, pickle, lifelines


Codes :
1. Analysis
- Survival analysis on the cohort of patients who underwent BCR.
The following R script performs Kaplan-Meier estimation of survival curves of the patients in the cohort across different risk groups. It also performs the Cox Regression on the cohort to identify clinical features predictive of BCR.
Latest codes : /survival_analysis_bcr/cox_reg_km_plotter_bcr.R

- Survival analysis on the cohort of patients with prostate biopsy pathology records demonstrating malignancy 
The following R script performs Kaplan-Meier estimation of survival curves of the patients in the cohort. It also performs the covariate adjustment and propensity re-weighting using Cox Proportional Hazard refression framework. 
Latest code : /survival_analysis_biopsy/ipw_survival.R
- Estimation of Heterogeneous group treatment effects
The following code performs propensity score prediction which is to be used for re-weighting hazard functions in Cox Regression. It also performs visualization of estimated survival curves across different risk groups. Finally, using Double Machine Learning framework (https://docs.doubleml.org/stable/index.html), it estimates average treatment effect as well as group average treatment effects across the risk groups.  
Latest code : /double_ml_causal_analysis/prostate_causal_inference.ipynb

2. Pre-processing
- Rule-based NLP algorithm for extracting clinical features from the unstructured notes. See https://docs.google.com/document/d/1aUFjvz8bumhCnSUDJ8w4CtKoNmGBKYTGeq0cMqoOOl0/edit and https://docs.google.com/document/d/1pRA2XcAxjbqbji8WHFfJD2bvwgo-EgryGLNRiq-gYT4/edit for more details on how the algorithm works on Prostate biopsy pathology report and Radical prostatectomy Pathology report, respectively.
  Latest codes : pre_processing/process_pat_v4.py (pathology report) and pre_processing/process_pat_biopsy.py (biopsy report)
- Data merging script. It merges the processed files from different sources of EHR data. 
  Latest code : pre_processing/data_merger.py
- Survival Data generator
  The following codes create survival data for the BCR survival analysis and causal inference analysis. Latest codes : pre_processing/create_cox_df.py, pre_processing/create_cox_df_causal.py 
  


3. Processed data
- Cannot be included due to the proprietary nature of the data (Need IRB approval)
