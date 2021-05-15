# Using Machine Learning to Uncover Prognosis of Patients with Localized Prostate Cancer

Abstract : 
Prostate cancer is one of the most commonly diagnosed cancer among males in the US. Although survival rates of patients diagnosed with localized disease are fairly high following Radical Prostatectomy (RP), about 10 - 15 \% patients experience biochemical recurrence (BCR). To gain clinical insights into the prognosis of patients with the disease, it is crucial to utilize Electronic Health Records (EHR) to build a data based clinical decision support system. To that end, we built the rule-based NLP algorithm which extracts clinical features in the unstructured pathology notes, obtained from Research Patient Data Registry, a centralized clinical data registry that gathers clinical information from various hospital systems. We validated our model's performance with a clinician and the model achieved overall 90\% accuracy on all abstracted features from 20 randomly selected pathology reports. Using those features along with other clinical information including age, ethnicity, and comorbidity, we performed survival analysis on the cohort of prostate cancer patients (n = 1,751) who underwent RP and discovered the extracted clinical features such as overall grade group and pathological tumor stages were significantly associated with the increased risk of BCR. Furthermore, we used the doubly-robust causal inference framework to estimate treatment effect of RP/Radiation on the survival of the cohort of patients with prostate biopsy pathology records demonstrating malignancy (n = 3,945). We estimates heterogeneous group treatment effects across different risk groups defined by age, overall disease stage, and comorbidities and found the strong evidence of survival benefit of the treatment across the groups. 


Codes :
1. Analysis
2. Pre-processing
3. Processed data
