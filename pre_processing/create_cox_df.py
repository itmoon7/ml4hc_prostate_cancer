import pandas as pd
import numpy as np
from tqdm import tqdm_notebook
import pickle

from lifelines.datasets import load_rossi
from lifelines import CoxPHFitter


def main(outcome_oi):
	# get RP procedure date
	with open('../data/empi_to_rp_date_dic.pkl', 'rb') as handle:
		empi_to_rp_date_dic = pickle.load(handle)
	# rp date min : 1992-8-5, rp date max : 2020-7-29
	rp_date_range = [np.arange(1992, 1998), np.arange(1998, 2004), np.arange(2004, 2009), np.arange(2009, 2015), np.arange(2015, 2021)]

	# get relevant data
	data_merged = pd.read_csv('../data/merged/df_merged_rp_positive_based_E.csv')
	data_merged.set_index('Unnamed: 0', inplace = True)
	data_merged.index.name = 'EMPI'

	data_merged_comorb = pd.read_csv('../data/merged/df_merged_rp_positive_based_C.csv')
	data_merged_comorb.set_index('Unnamed: 0', inplace = True)
	data_merged_comorb.index.name = 'EMPI'

	empis_oi = set(data_merged.index) & set(data_merged_comorb.index)
	comorbs_oi = ['ami_agg', 'chf_agg', 'pvd_agg', 'cevd_agg', 'copd_agg', 'rheumd_agg', 'pud_agg', 'mld_agg', 'diab_agg', 'diabwc_agg', 'hp_agg', 'rend_agg', 'metacanc_agg', 'aids_agg']
	data_merged = pd.concat([data_merged.loc[empis_oi], data_merged_comorb.loc[empis_oi][comorbs_oi]], axis = 1)
	# process pt stage
	data_merged['pt1'] = 0; data_merged['pt1a'] = 0; data_merged['pt1b'] = 0; data_merged['pt1c'] = 0
	data_merged['pt2'] = 0; data_merged['pt2a'] = 0; data_merged['pt2b'] = 0; data_merged['pt2c'] = 0
	data_merged['pt3'] = 0; data_merged['pt3a'] = 0; data_merged['pt3b'] = 0; data_merged['pt3c'] = 0
	data_merged['pt4'] = 0

	for empi, pt_stage in zip(data_merged.index, data_merged.pT_stage_combined.values):
		data_merged.at[empi, pt_stage] = 1
	#     if pt_stage[:-1] in data_merged.columns:
	#         data_merged.at[empi, pt_stage[:-1]] = 1
	# data_merged.drop(columns = ['pT_stage_combined'], inplace = True)

	# get rp procedure date
	rp_date_list = []
	for empi in data_merged.index:
		rp_date = empi_to_rp_date_dic[empi]
		for yr_idx, yr_range in enumerate(rp_date_range):
			if rp_date.year in yr_range:
				rp_date_list.append(yr_idx)
				break
	data_merged['rp_date'] = rp_date_list

	# load outcome : 
	df_outcome = pd.read_csv('../data/df_outcome_rp_positive_final.csv')
	df_outcome.set_index('EMPI', inplace = True)

	# outcome_oi = 'death' # 'death', 'bcr'
	empis_oi = set(df_outcome.index) & set(data_merged.index)
	df_outcome_oi = df_outcome.loc[empis_oi]
	df_merged_data_oi = data_merged.loc[empis_oi]

	print('Feature stats in the cohort of interest : ')
	print(df_merged_data_oi.sum())

	drop_pt_stage = ['pt1', 'pt1a', 'pt1b', 'pt1c', 'pt3c', 'pt4', 'pt3']# 'Unknown/other', 'pt3', 'Asian']
	drop_other = ['Unknown/other', 'Asian']
	drop_cols = drop_pt_stage + drop_other
	df_merged_data_oi.drop(columns = drop_cols, inplace = True)
	# df_merged_data_oi = df_merged_data_oi.loc[~df_merged_data_oi.isin(drop_pt_stage)]


	print('New featur stats in the cohort of interest : ')
	print(df_merged_data_oi.sum())

	if outcome_oi == 'bcr':
		outcome_cols_oi = ['bcr_ind', 'time_to_bcr_in_month']
		df_outcome_oi = df_outcome_oi[outcome_cols_oi]
	else:
		outcome_cols_oi = ['death_ind', 'time_to_death_in_month']
		df_outcome_oi = df_outcome_oi[outcome_cols_oi]
	df_cox_data = pd.concat([df_merged_data_oi, df_outcome_oi], axis = 1)
	print('\n')
	print('Exporting data...')
	if outcome_oi == 'bcr':
		df_cox_data.to_csv('../data/df_cox_data_bcr.csv')
	else:
		df_cox_data.to_csv('../data/df_cox_data_death.csv')
	print('\n')
	run_cox = True
	if run_cox:
		cph = CoxPHFitter(l1_ratio=1, penalizer = 0.01)
		if outcome_oi == 'death':
			cph.fit(df_cox_data.drop(columns = ['pT_stage_combined']), 'time_to_death_in_month', 'death_ind')
			cph.print_summary()
			print(cph.summary)
			cph.summary.round(3).to_csv('cox_result_death.csv')
		else:
			cph.fit(df_cox_data.drop(columns = ['pT_stage_combined']), 'time_to_bcr_in_month', 'bcr_ind')
			cph.print_summary()
			print(cph.summary)
			cph.summary.round(3).to_csv('cox_result_bcr.csv')
		breakpoint()
		"""
		drop_cols = ['max_psa', 'min_psa', 'mean_psa', 'ami_agg', 'chf_agg', 'pvd_agg', 'cevd_agg', 'copd_agg', 'rheumd_agg', 'pud_agg', 'mld_agg', 'diabwc_agg', 'hp_agg', 'rend_agg', 'metacanc_agg', 'aids_agg']
		"""
	return

if __name__ == "__main__":
	outcome_oi = input('Choose the outcome of interest (bcr or death) : ')
	if outcome_oi not in ['bcr', 'death']:
		raise KeyError('Only supports bcr or death atm')
	main(outcome_oi)