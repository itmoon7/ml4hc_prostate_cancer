import pandas as pd
import numpy as np
from tqdm import tqdm_notebook
import pickle

from lifelines.datasets import load_rossi
from lifelines import CoxPHFitter


def main(filter_gap_days, rp_vs_radiation):
	# get RP procedure date
	with open('../data/empi_to_rp_date_dic.pkl', 'rb') as handle:
		empi_to_rp_date_dic = pickle.load(handle)

	# if rp_vs_radiation:
	# get biopsy date
	with open('../data/empi_to_date_oi_dic.pkl', 'rb') as handle:
		empi_to_date_oi_dic = pickle.load(handle)
	# rp date min : 1992-8-5, rp date max : 2020-7-29
	# rp_date_range = [np.arange(1992, 1998), np.arange(1998, 2004), np.arange(2004, 2009), np.arange(2009, 2015), np.arange(2015, 2021)]

	# pre-process the data using radiation only patients / multiple RP patients 
	df_first_date_rads = pd.read_csv('../data/processed_data/first_date_rads.csv')
	df_first_date_rads.set_index('EMPI', inplace = True)
	for empi, pr_date in zip(df_first_date_rads.index, df_first_date_rads.prdate_parsed.values):
		if empi in empi_to_date_oi_dic.keys():
			df_first_date_rads.at[empi, 'pr_date_minus_biopsy_date'] = (pd.to_datetime(pr_date) - empi_to_date_oi_dic[empi]).days
	df_first_date_rads_filtered = df_first_date_rads.loc[df_first_date_rads.pr_date_minus_biopsy_date > 0]
	df_first_date_rads_filtered = df_first_date_rads_filtered.loc[df_first_date_rads_filtered.pr_date_minus_biopsy_date <= filter_gap_days]
	radiation_empis = set(df_first_date_rads_filtered.index)

	df_multirp = pd.read_csv('../data/processed_data/multirp.csv')
	df_multirp.set_index('EMPI', inplace = True)
	multirp_empis = set(df_multirp.index)
	# breakpoint()

	# get relevant data
	df_merged_comorb = pd.read_csv('../data/merged/df_merged_biopsy_based_C.csv')
	df_merged_comorb.set_index(df_merged_comorb.columns[0], inplace = True)
	df_merged_comorb.index.name = 'EMPI'

	df_merged_psa_prior = pd.read_csv('../data/merged/df_merged_biopsy_based_E.csv')
	df_merged_psa_prior.set_index(df_merged_psa_prior.columns[0], inplace = True)
	df_merged_psa_prior.index.name = 'EMPI'

	df_outcome_oi = pd.read_csv('../data/df_outcome_final_biopsy_based_clean.csv')
	df_outcome_oi.set_index('EMPI', inplace = True)

	# filter_gap_days = 60
	df_outcome_oi_filtered_rp_postive = df_outcome_oi.loc[df_outcome_oi.rp_date_minus_biopsy_date_in_days <= filter_gap_days]
	df_outcome_oi_rp_negative = df_outcome_oi.loc[df_outcome_oi.rp_date.isnull()]
	if rp_vs_radiation:
		# exclude both RP and Radiation
		both_rp_radiation_empis = set(df_outcome_oi_filtered_rp_postive.index) & radiation_empis
		df_outcome_oi_filtered_rp_postive = df_outcome_oi_filtered_rp_postive.loc[set(df_outcome_oi_filtered_rp_postive.index) - both_rp_radiation_empis]
		df_outcome_oi_rp_negative = df_outcome_oi_rp_negative.loc[set(df_outcome_oi_rp_negative.index) & radiation_empis - multirp_empis]
		df_outcome_oi_filtered_total = pd.concat([df_outcome_oi_filtered_rp_postive, df_outcome_oi_rp_negative])
	else: # treated (radiation and rp) vs AS
		df_outcome_oi_treated_postive = pd.concat([df_outcome_oi_filtered_rp_postive, df_outcome_oi_rp_negative.loc[set(df_outcome_oi_rp_negative.index) & radiation_empis]])
		# df_outcome_oi_filtered_rp_postive
		df_outcome_oi_treated_negative = df_outcome_oi_rp_negative.loc[set(df_outcome_oi_rp_negative.index) - radiation_empis]
		df_outcome_oi_filtered_total = pd.concat([df_outcome_oi_treated_postive, df_outcome_oi_treated_negative])

	# filter out negative time to death
	df_outcome_oi_filtered_final = df_outcome_oi_filtered_total.loc[df_outcome_oi_filtered_total.time_to_death_in_month > 0]
	# df_outcome_oi_filtered_final.set_index('EMPI', inplace = True)

	# merge comorbidity data and psa prior data
	empi_common = set(df_merged_comorb.index) & set(df_merged_psa_prior.index)
	df_merged_psa_comorb = pd.concat([df_merged_comorb.loc[empi_common], df_merged_psa_prior.loc[empi_common].psa_prior_to_rp], axis = 1)
	df_merged_psa_comorb.drop(columns = ['wscore_agg'], inplace = True)

	# merge outcome and feature dfs
	empi_common_outcome = set(df_outcome_oi_filtered_final.index) & set(df_merged_psa_comorb.index)
	outcome_cols_oi = ['death_ind', 'time_to_death_in_month']
	df_cox = pd.concat([df_merged_psa_comorb.loc[empi_common_outcome], df_outcome_oi_filtered_final.loc[empi_common_outcome][outcome_cols_oi]], axis = 1)

	if not rp_vs_radiation:
		df_cox.loc[df_cox.index.isin(df_outcome_oi_treated_postive.index), 'rp_indicator'] = 1
		df_cox.rename(columns = {'rp_indicator':'treated'}, inplace = True)
	# get biopsy date range
	# biopsy_date_list = []
	# for empi in data_merged.index:
	# 	biopsy_date = empi_to_rp_date_dic[empi]
	# 	for yr_idx, yr_range in enumerate(rp_date_range):
	# 		if biopsy_date.year in yr_range:
	# 			biopsy_date_list.append(yr_idx)
	# 			break
	# data_merged['rp_date'] = biopsy_date_list



	# drop uninformative and extreme minority features
	drop_cols = ['benign', 'Unknown/other'] # pre-filtering
	# remove metastatic cancer patients
	df_cox = df_cox.loc[df_cox.metacanc_agg == 0]
	# drop_cols_non_informative = ['cevd_agg', 'rheumd_agg', 'pud_agg', 'mld_agg', 'diabwc_agg', 'aids_agg', 'metacanc_agg']
	if rp_vs_radiation:
		drop_cols_non_informative = ['metacanc_agg', 'rheumd_agg', 'copd_agg', 'mld_agg', 'diabwc_agg', 'hp_agg', 'rend_agg', 'aids_agg']
	else:
		drop_cols_non_informative = ['metacanc_agg', 'cevd_agg', 'rheumd_agg', 'mld_agg', 'diabwc_agg', 'hp_agg', 'aids_agg']#, 'rheumd_agg', 'copd_agg', 'mld_agg', 'diabwc_agg', 'hp_agg', 'rend_agg', 'aids_agg']
	df_cox_final = df_cox.drop(columns = drop_cols + drop_cols_non_informative)
	# standardize age and auxiiliary_mci_score,psa_prior_to_rp
	standardize_cols = ['Age at RP', 'auxiiliary_mci_score', 'psa_prior_to_rp']
	for col in standardize_cols:
		if col == 'psa_prior_to_rp':
			breakpoint()
		df_cox_final[col] = (df_cox_final[col].values - np.mean(df_cox_final[col]))/np.std(df_cox_final[col].values)

	
	# df_cox_final = df_cox_final.loc[df_cox_final.overall_grade_merged == 1]
	# df_cox_final.drop(columns = ['overall_grade_merged'], inplace = True)
	# df_cox_final = df_cox_final.loc[df_cox_final.overall_grade_merged > 1]
	# df_cox_final.drop(columns = ['overall_grade_merged'], inplace = True)

	print('Final cox df stats : ')
	print(df_cox_final.sum())

	cph = CoxPHFitter(penalizer = 0.00, l1_ratio = 0)
	cph.fit(df_cox_final, 'time_to_death_in_month', 'death_ind', show_progress = False, step_size = 0.1)
	cph.print_summary()

	breakpoint()

	# rename colums for downstraem task compatibility
	df_cox_final.rename(columns = {'death_ind':'death', 'rp_indicator':'rp', 'time_to_death_in_month':'survtime'}, inplace = True)
	if rp_vs_radiation:
		df_cox_final.to_csv('../data/df_cox_data_death_causal_inference_rp_vs_radiation.csv')
	else:
		df_cox_final.to_csv('../data/df_cox_data_death_causal_inference.csv')
	breakpoint()
	return

if __name__ == "__main__":
	"""
	In one causal analysis, we want to estimate CATE among patients undergoing RP vs. radiation. The most clinically relevant question is in those with grade group 2+
	In another causal analysis, we want estimate CATE among patients undergoing RP vs. radiation vs. AS. The most clinically relevant question is in those with grade group 1.
	In another causal analysis (there are many possible question), we want to estimate CATE among patients with BCR after RP comparing RP + salvage radiation vs. RP + other treatment/no treatment
	"""
	filter_gap_days = float(input('Choose gap days between RP date and biopsy date : '))
	rp_vs_radiation = input('RP vs. radiation? (Y/N) : ') == 'Y'
	main(filter_gap_days, rp_vs_radiation)