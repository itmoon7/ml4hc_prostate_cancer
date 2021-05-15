import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import re

"""
Intae Moon
Date : April 3, 2021
"""
# simple function to detect negation

def featurize_biopsy_df(df_biopsy):
	# get the following features 
	# 1. mci auxiiliary feature
	df_biopsy_wo_duplicate = df_biopsy.drop_duplicates(subset = 'EMPI', keep = False)
	df_biopsy_wo_duplicate.set_index('EMPI', inplace = True)
	empis_oi = np.unique(df_biopsy_wo_duplicate.index)
	counter = 0
	for empi in empis_oi:
		num_pos_cores_oi = df_biopsy_wo_duplicate.at[empi, 'num_pos_cores']
		overall_grade_group_oi = df_biopsy_wo_duplicate.at[empi, 'overall_grade_group']
		mci_oi = df_biopsy_wo_duplicate.at[empi, 'max_core_involve']
		if len(num_pos_cores_oi) == len(overall_grade_group_oi) == len(mci_oi):
			weighted_val = 0
			for num_pos_core, overall_grade, mci in zip(num_pos_cores_oi, overall_grade_group_oi, mci_oi):
				if len(mci) == num_pos_core: # mci = [a,b,c], num_pos_core = 3, result = a + b + c
					for mci_sub in mci:
						weighted_val += 1 * overall_grade * mci_sub/100 # 
				else: # mci = [a,b], num_pos_core = 3, result = a + b + b or mci = [a,b,c,d], num_pos_core = 3, result = a + b + c
					# if len(mci) == 1:
					# 	weighted_val += 1 * overall_grade * mci_sub/100/num_pos_core # mci = [a], 
					# else:
					for i in range(num_pos_core):
						try:
							weighted_val += 1 * overall_grade * mci[i]/100
						except:
							weighted_val += 1 * overall_grade * mci[-1]/100
			df_biopsy_wo_duplicate.at[empi, 'auxiiliary_mci_score'] = weighted_val
		else:
			counter += 1
	df_biopsy_wo_duplicate.dropna(subset = ['auxiiliary_mci_score'], inplace = True)
	return df_biopsy_wo_duplicate

def get_tumor_stage_info_via_context(df_pathology_with_biopsy):
	tumor_stage_list = []
	hard_coded_pt2_stage = ["Extraprostatic extension of tumor is not identified", "confined to the prostate", "Extraprostatic extension of tumor is not identified", "confined to the gland", "confined to prostate", "confined within the prostatic capsule", "involves the capsule but is confined to the prostate", "invades but does not transgress the capsule", "No capsular penetration demonstrated", "NOT INFILTRATING PERIPROSTATIC ADIPOSE TISSUE", "tumor does not transgress the prostatic capsule", "apparently localized", "invades into but not through the prostatic capsule", "extends to, but not through the prostatic capsule.", "the periprostatic and soft tissue resection margins are free of carcinoma", "no capsular invasion is present.", "extracapsular extension is not identified", "confined within the prostatic capsule", "TUMOR IS NOT SEEN OUTSIDE THE CAPSULE", "does not extend into periprostatic fat", "is confirmed to the prostate", "Extraprostatic extension is not present", "tumor does not transgress the prostatic capsule", "All inked margins and capsule are free of tumor", "no tumor involvement of the capsule", "Extraprostatic extension is not identified", "No extraprostatic extension of tumor"]
	hard_coded_pt2_stage = [val.lower() for val in hard_coded_pt2_stage]

	hard_coded_pt3a_stage = ["extends into the prostatic", "extends through the capsule into periprostatic", "extends slightly into periprostatic tissue", "LIMITED EXTRAPROSTATIC EXTENSION", "tumor extends into extraprostatic soft tissue", "tumor focally penetrates through the prostate capsule", "focally into the periprostatic fat"]
	hard_coded_pt3a_stage = [val.lower() for val in hard_coded_pt3a_stage]

	hard_coded_pt3b_stage = ["extends into extra prostatic soft tissues and the left seminal vesicle", "and involves the left seminal vesicle"]
	hard_coded_pt3b_stage = [val.lower() for val in hard_coded_pt3b_stage]
	
	for text_oi, report_number in tqdm(zip(df_pathology_with_biopsy.Report_Text.values, df_pathology_with_biopsy.Report_Number.values), total = len(df_pathology_with_biopsy), desc = 'extracing pTstage via context...'):
		# lower case all the words in each report
		text_oi_list = text_oi.lower().split(' ')
		tumor_stage_ind = 0; tumor_stage_ind_pt3a = 0; tumor_stage_ind_pt2 = 0; tumor_stage_ind_pt3b = 0

		# check for hard-coded rules (pT stage)
		for val_harcoded_pt2 in hard_coded_pt2_stage:
			if val_harcoded_pt2 in text_oi.lower():
				tumor_stage_ind = 1
				tumor_stage_list.append('pt2')
				break

		for val_harcoded_pt3a in hard_coded_pt3a_stage:
			if val_harcoded_pt3a in text_oi.lower():
				if tumor_stage_ind:
					print('(pT stage) presence of conflicting sentences in Report : ', report_number)
					print(val_harcoded_pt2)
					print(val_harcoded_pt3a)
					print('\n')
					tumor_stage_list.pop()
					tumor_stage_list.append('pt3a')
					# breakpoint()
				else:
					tumor_stage_ind = 1
					tumor_stage_list.append('pt3a')
				break

		for val_harcoded_pt3b in hard_coded_pt3b_stage:
			if val_harcoded_pt3b in text_oi.lower():
				if tumor_stage_ind:
					print('(pT stage) presence of conflicting sentences in Report : ', report_number)
					print(val_harcoded_pt3a)
					print(val_harcoded_pt3b)
					print('\n')
					tumor_stage_list.pop()
					tumor_stage_list.append('pt3b')
					# breakpoint()
				else:
					tumor_stage_ind = 1
					tumor_stage_list.append('pt3b')
				break


		for idx, val in enumerate(text_oi_list):
			# get tumor stage info
			if not tumor_stage_ind:
				if 'confine' in val: # check for "confined to the prostate"
					# check for negation :
					if 'not' not in text_oi_list[idx -1]:
						context_to_look = text_oi_list[idx : idx + 5]
						found_to = 0
		#                 for query_word in context_to_look:
						for context_word in context_to_look:
							if 'to' in context_word:
								found_to = 1
							if 'prostate' in context_word and found_to:
								tumor_stage_ind = 1
								tumor_stage_ind_pt2 = 1
								tumor_stage_list.append('pt2')
		#                         print('(Hit) confined to prostate')
		#                         print(context_to_look)
								break
				elif 'extend' in val:
					idx_extend = idx
					context_to_look = text_oi_list[idx : idx + 15]
					found_to = 0
					for idx, context_word in enumerate(context_to_look):
						if 'to' in context_word: # extends "to" xxx : 
							found_to = 1
						if 'soft' in context_word and found_to:
							if idx < len(context_to_look) - 1:
								if 'tissue' in context_to_look[idx+1]:
									tumor_stage_ind = 1
									tumor_stage_ind_pt3a = 1
									tumor_stage_list.append('pt3a')
									break

					found_to = 0
					for idx, context_word in enumerate(context_to_look):
						if 'to' in context_word: # extends "to" xxx : 
							found_to = 1
						if 'prostatic' in context_word and found_to:
							if idx < len(context_to_look) - 1:
								if 'capsule' in context_to_look[idx+1]:
									if not tumor_stage_ind_pt3a:
										tumor_stage_ind = 1
										tumor_stage_ind_pt3a = 1
										tumor_stage_list.append('pt3a')
	#                                 print('(Hit) extend to prostatic capsule')
	#                                 print(context_to_look)
									break
					if 'negative' not in context_to_look and 'no' not in context_to_look and 'not' not in context_to_look: #check for negation of seminal vesicle
						found_to = 0
						for idx, context_word in enumerate(context_to_look):
							if 'to' in context_word: # extends "to" xxx : 
								found_to = 1
							if 'seminal' in context_word and found_to:
								if idx < len(context_to_look) - 1:
									if 'vesicle' in context_to_look[idx+1] or 'vesicles' in context_to_look[idx+1]:
										if tumor_stage_ind_pt3a == 1:
											# print('duplicate assignment!')
											# print('fixing it...')
											tumor_stage_list.pop()
											tumor_stage_list.append('pt3b')
											tumor_stage_ind = 1
											tumor_stage_ind_pt3b = 1
										else:
											tumor_stage_list.append('pt3b')
											tumor_stage_ind = 1
											tumor_stage_ind_pt3b = 1
										break
					# if report_number == 'BS13K52639':
					# 	breakpoint()
			if tumor_stage_ind_pt3a == 1: # fix pt3a case if you find seminal later in the report
				if 'seminal' in val:
					if idx - idx_extend > 35: # more than 35 words betweeen extend and seminal, which means seminal info is irrelevant
						break
					context_to_look = text_oi_list[idx - 7 : idx + 7]
					if 'negative' not in context_to_look and 'no' not in context_to_look and 'not' not in context_to_look:
						tumor_stage_list.pop()
						tumor_stage_list.append('pt3b')
					else:
						tumor_stage_list.pop()
						tumor_stage_list.append('pt3a')
					break # break it here since you only want to take a look at seminal context right after extend
		# update tumor stage
		if not tumor_stage_ind:
			tumor_stage_list.append(None)
	df_pathology_with_biopsy['pT_stage_context'] = tumor_stage_list
	return df_pathology_with_biopsy


def func_det_negation(text_detect_neg):
	det_negation = 0
	for idx, val in enumerate(text_detect_neg):
		if 'no' == val.lower():
			bool_list = []
			for val_ in text_detect_neg[idx:]:
				if '.' not in val_:
					bool_list.append(False)
				else:
					bool_list.append(True)
			if not any(bool_list): # if there is no dot between 'no' --- 'query word', then we consider it as a negation
				det_negation = 1
	return det_negation
	
def main():
	load_prcossed = input('Load the processed data? (Y/N) : ') == 'Y'
	if load_prcossed:
		df_pathology_with_biopsy = pd.read_csv('df_pathology_with_biopsy_params_extracted.tsv', sep = '\t')
		df_pathology_with_biopsy.set_index(df_pathology_with_biopsy.columns[0], inplace = True)
	else:
		# load the processed pathology report
		df_pathology = pd.read_csv('../data/Pat_total.csv', sep="|", low_memory=False)
		df_pathology = df_pathology.loc[df_pathology.Report_Description == 'Surgical Pathology']

		# get all the pathology reports with radical prostatectomy included
		idx_oi = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'needle' in report.lower() or 'core' in report.lower() or 'biops' in report.lower() and 'prost' in report.lower()]
		idx_to_exclude = [idx for idx, report in enumerate(df_pathology.Report_Text.values) if 'prostatectomy' in report.lower()]# or 'nephrectomy' in report.lower() or 'cystectomy' in report.lower() or 'cystoprostatectomy' in report.lower() or 'ureter' in report.lower()] # or 'prostate, radical resection' in report.lower()]

		df_pathology_with_biopsy = df_pathology.iloc[np.sort(list(set(idx_oi) - set(idx_to_exclude)))].copy()
		# loading singlerp.csv 
		single_rp_patients = pd.read_csv('../data/singlerp.csv')
		# df_pathology_with_biopsy = df_pathology_with_biopsy.loc[df_pathology_with_biopsy.EMPI.isin(single_rp_patients.EMPI.values)]

		directional_words = ['right', 'left']
		primary_grade_list = []; secondary_grade_list = []; overall_grade_list = []; overall_grade_merged_list = []; overall_gs_list = [];
		radical_pros = []; tumor_stage_list = []; margin_list = []; benign_list = []
		num_pos_cores_list = []; num_pos_cores_sum_list = []; num_total_core_list = []; num_total_core_sum_list = []; max_core_involve_list = [];
		small_cell_carc = []; neuroendocrine_carc = []; adenocarcinoma = []
		stage_count = 0; margin_counter = 0; margin_counter_hard_coded = 0; pt_stage_counter_hard_coded = 0

		# for num core involvement 		
		num_words_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
		desc_num_words_list = ['one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten']
		desc_num_to_num_dic = {desc_num:num for num, desc_num in zip(num_words_list, desc_num_words_list)}
		for text_oi, report_number in tqdm(zip(df_pathology_with_biopsy.Report_Text.values, df_pathology_with_biopsy.Report_Number.values), total = len(df_pathology_with_biopsy), desc = 'extracing info...'):
			# create local grade lists
			primary_grade_list_local = []; secondary_grade_list_local = []; overall_grade_list_local = []; overall_gs_list_local = []
			num_pos_core_list_local = []; num_total_core_list_local = []; max_core_involve_list_local = [];
			num_pos_cores = 0; radical_pros_ind = 0; gleason_reported = 0; 
		#   all patient should have received rp
			radical_pros_ind = 1
			radical_pros.append(1)
			
			# sentence_to_check = "INVOLVEMENT OF MARGIN"
			text_oi = text_oi.lower()

			# benign tumor? 
			if 'gleason score' in text_oi or 'gleason grade' in text_oi:
				benign_list.append(0)
			else:
				if ('no tumor seen' in text_oi or 'no carcinoma' in text_oi):
					benign_list.append(1)
				else:
					benign_list.append(None)

			text_oi_list = text_oi.split(' ')
			found_directional_word = 0; margin_ind = 0
			small_cell_carc_ind = 0; neuroendocrine_carc_ind = 0; adenocarcinoma_ind = 0; tumor_stage_ind = 0
			
			for idx, val in enumerate(text_oi_list):
				if any(dir_val in val for dir_val in directional_words):
					found_directional_word = 1
				
				# counts number of positive cores
				# if 'gleason' in val and found_directional_word:
				# 	if 'm.d.,' not in text_oi_list[idx - 2 : idx + 2]: # filter out any MDs with last name Gleason
				# 		num_pos_cores += 1
					
				# histology : 
				# small cell carcinoma
				if 'small' in val and 'cell' in text_oi_list[idx+1]:
					if 'carcinoma' in [txt for txt in text_oi_list[idx + 1:idx + 5]]:
						text_detect_neg = text_oi_list[idx - 3:idx]
						det_negation = func_det_negation(text_detect_neg)
						if not det_negation:
							small_cell_carc_ind = 1
		#                 else:
		#                     print('small-cell', text_oi_list[idx - 5:idx + 5])
		#                     breakpoint()
				# neuroendocrine carcinoma
				if 'neuroendocrine' in val and 'carcinoma' in text_oi_list[idx + 1]:
					text_detect_neg = text_oi_list[idx - 3:idx]
					det_negation = func_det_negation(text_detect_neg)
					if not det_negation:
						neuroendocrine_carc_ind = 1
		#             else:
		#                 print('neuroendocrine', text_oi_list[idx - 5:idx + 5])
		#                 breakpoint()
				# adenocarcinoma
				if 'adenocarcinoma' in val:
					text_detect_neg = text_oi_list[idx - 3:idx]
					det_negation = func_det_negation(text_detect_neg)
					if not det_negation:
						adenocarcinoma_ind = 1
		#             else:
		#                 print('adenocarcinoma', text_oi_list[idx - 5:idx + 5])
		#                 breakpoint()
				# if gleason score comes before any biopsy information, indicated by the directional words like "right" base, "left" base, and etc.
				# we report overall grade, each biopsy grade, and number of positive cores / total cores 
				primary_grade = None; secondary_grade = None; num_pos_core = None; num_total_core = None; max_core_involve = [];
				# if 'gleason' in val and report_number == 'S13-36632':
				# 	print(text_oi_list[idx:idx + 15]) # heuristically choose 7. scores should be included in the next 7 wor)
				# 	print('Missed?')
				# 	breakpoint()
				if 'gleason' in val and found_directional_word and 'm.d.,' not in text_oi_list[idx - 2 : idx + 2]:# not found_directional_word: 
					# get primary and secondary grades
					found_primary_grade = 0; found_secondary_grade = 0; found_plus_sign = 0; found_one_dash = 0
					raw_gleason_text = ''.join(text_oi_list[idx:idx + 15]) # heuristically choose 7. scores should be included in the next 7 words
					
					# check number of positive and total cores
					raw_gleason_cores_contex_joined = ' '.join(text_oi_list[idx:idx + 25])
					raw_gleason_cores_context  = text_oi_list[idx:idx + 25]
					if 'of' in raw_gleason_cores_context:
						for idx_core_word, core_word in enumerate(raw_gleason_cores_context):
							if core_word == 'of':# and :
								context_core_oi = raw_gleason_cores_context[idx_core_word-1:idx_core_word+2]
								# if successfully extracted three words
								if len(context_core_oi) == 3:
									pos_core_word_oi = re.sub('[()]', '', context_core_oi[0])
									total_core_word_oi = re.sub('[()]', '', context_core_oi[2])								
									# involving num of num
									if pos_core_word_oi.isdigit() and total_core_word_oi.isdigit():
										num_pos_core = int(pos_core_word_oi)
										num_total_core = int(total_core_word_oi)
									# involving desc_num of desc_num
									elif pos_core_word_oi in desc_num_words_list and total_core_word_oi in desc_num_words_list:
										num_pos_core = desc_num_to_num_dic[pos_core_word_oi]
										num_total_core = desc_num_to_num_dic[total_core_word_oi]
									# involving desc_num (num) of desc_num
									elif pos_core_word_oi.isdigit() and total_core_word_oi in desc_num_words_list:
										num_pos_core = int(pos_core_word_oi)
										num_total_core = desc_num_to_num_dic[total_core_word_oi]
									# or vice versa
									elif pos_core_word_oi in desc_num_words_list and total_core_word_oi.isdigit():
										num_pos_core = desc_num_to_num_dic[pos_core_word_oi]
										num_total_core = int(total_core_word_oi)
									if type(num_total_core) == int:
										break
										# print(raw_gleason_cores_context)
										# print(context_core_oi)
										# print(num_pos_core)
										# print(num_total_core)
										# breakpoint()

					# get MCI metrics
					if '%' in raw_gleason_cores_contex_joined and num_pos_core is not None:
						for idx_mci, char_mci in enumerate(raw_gleason_cores_contex_joined):
							if char_mci == '%' and len(max_core_involve) < num_pos_core:
							# percent_index = raw_gleason_cores_contex_joined.find('%')
							# print()
								percent_word = raw_gleason_cores_contex_joined[idx_mci-3:idx_mci+1]
								# in this scenario 5-10% only get 10%
								percent_word = percent_word.split('-')[-1]
								# take only numeric parts
								percent_numeric = re.sub('\D', '', percent_word)
								if percent_numeric.isdigit():
									if percent_numeric == '00':
										max_core_involve.append(100)
									else:
										max_core_involve.append(int(percent_numeric))
						# if len(max_core_involve) > 1:
						# 	print(raw_gleason_cores_context)
						# 	print(num_pos_core, num_total_core)
						# 	print(max_core_involve)
						# 	breakpoint()
						# pass

					# remove text in () or [] : e.g. Score 3 (60%) + 4 (40%) = 7/10 -> Score 3 + 4 = 7/10
					raw_gleason_text = re.sub("[\(\[].*?[\)\]]", "", raw_gleason_text)
					gleason_reported = 1 # this ensures we only look at the first appearance of gleason score

					# check if there's only one dash : to capture gleason 3/5 where primary 3 and secondary 3
					slash_counter = 0; dash_counter = 0; plus_sign_counter = 0; only_one_slash = 0
					raw_gleason_text_slash = ''.join(text_oi_list[idx:idx + 15]) # capture a little longer context
					for char_dash in raw_gleason_text_slash:
						if '/' in char_dash:
							slash_counter += 1
						if '+' in char_dash:
							plus_sign_counter += 1
						if '-' in char_dash:
							dash_counter += 1

					if slash_counter == 1 and plus_sign_counter == 0 and dash_counter == 0:
						only_one_slash = 1
						# print(raw_gleason_text_slash)
						

					for idx_char, char in enumerate(raw_gleason_text):
						# overall grade must be represented as primary + secondary in the report
						if char == '+':
							raw_gleason_equation = raw_gleason_text[idx_char-10:idx_char+10]
							for sub_char in raw_gleason_equation:
								# checking + sign again for a few corner cases 
								if sub_char == '+':
									found_plus_sign = 1
								if sub_char.isdigit():
									if not found_plus_sign:
										primary_grade_to_be = int(sub_char)
										found_primary_grade = 1
									elif found_primary_grade and found_plus_sign and not found_secondary_grade:
										primary_grade = primary_grade_to_be
										secondary_grade = int(sub_char)
										found_secondary_grade = 1
										break
							break
						elif char == '/' and not found_one_dash:
							# this case is to capture Gleason 3-4/5 
							raw_gleason_grade_dash = raw_gleason_text[idx_char-3:idx_char+2]
							if '-' in raw_gleason_grade_dash:
								raw_gleason_grade_dash_split = raw_gleason_grade_dash.split('/')
								# breakpoint()
								if raw_gleason_grade_dash_split[1] == '5' and raw_gleason_grade_dash_split[0][0].isdigit() and raw_gleason_grade_dash_split[0][2].isdigit():
									primary_grade = int(raw_gleason_grade_dash_split[0][0])
									secondary_grade = int(raw_gleason_grade_dash_split[0][2])
									# print(raw_gleason_grade_dash)
									break
							# this case is to capture Gleason 6/10 
							raw_gleason_grade_dash_and_ten = raw_gleason_text[idx_char-1:idx_char+3]
							raw_gleason_grade_dash_and_ten_split = raw_gleason_grade_dash_and_ten.split('/')
							if raw_gleason_grade_dash_and_ten_split[0].isdigit() and raw_gleason_grade_dash_and_ten_split[1] == '10':
								# then primary = secondary = 3 . This is not too important 
								primary_grade = int(raw_gleason_grade_dash_and_ten_split[0])/2
								secondary_grade = primary_grade
								# print(raw_gleason_text)
								# breakpoint()
								break

							# if it is the form of x/5 : gleason grade 3/5 and 2/5 confined to
							raw_gleason_grade = raw_gleason_text[idx_char-1:idx_char+2]
							raw_gleason_grade_split = raw_gleason_grade.split('/')
							if len(raw_gleason_grade_split) == 2 and raw_gleason_grade_split[0].isdigit() and raw_gleason_grade_split[1] == '5':
								found_one_dash = 1
								primary_grade_to_be = int(raw_gleason_grade_split[0])
								# print(text_oi_list[idx:idx + 7])
								# print(primary_grade_to_be)
								# breakpoint()
								continue # thie ensures the next if statement is true automatically 
    					# encounter second dash : most likely we have something like gleason grade 3/5 and 2/5 confined to...
						if char == '/' and found_one_dash:
							raw_gleason_grade = raw_gleason_text[idx_char-1:idx_char+2]
							raw_gleason_grade_split = raw_gleason_grade.split('/')
							# if it is the form of x/5
							if len(raw_gleason_grade_split) == 2 and raw_gleason_grade_split[0].isdigit() and raw_gleason_grade_split[1] == '5':
								primary_grade = primary_grade_to_be
								secondary_grade = int(raw_gleason_grade_split[0])
								# print(raw_gleason_equation)
								# print(text_oi_list[idx:idx + 7])
								# print('primary_grade : ', primary_grade)
								# print('secondary_grade : ', secondary_grade)
								# print('\n')
								break
						# check if there's only one dash : to capture gleason 3/5 where primary 3 and secondary 3
						elif found_one_dash and only_one_slash:
							primary_grade = primary_grade_to_be
							secondary_grade = primary_grade_to_be
							# print(raw_gleason_grade)
							# print(raw_gleason_text_slash)
							# print(text_oi_list[idx-10:idx + 25])
							# breakpoint()
							break

					# check
					# if primary_grade == 1:
					# 	print('Primary grade')
					# 	print(primary_grade)
					# 	print(raw_gleason_cores_context)
					# 	breakpoint()
					# if secondary_grade == 1:
					# 	print('Secondary grade')
					# 	print(secondary_grade)
					# 	print(raw_gleason_cores_context)
					# 	breakpoint()
					# if 595 in max_core_involve:
					# 	print('MCI')
					# 	print(max_core_involve)
					# 	print(raw_gleason_cores_context)
					# 	breakpoint()	
					# if 890 in max_core_involve:
					# 	print('MCI')
					# 	print(max_core_involve)
					# 	print(raw_gleason_cores_context)
					# 	breakpoint()						

				# https://www.pcf.org/about-prostate-cancer/diagnosis-staging-prostate-cancer/gleason-score-isup-grade/
				if primary_grade is not None and secondary_grade is not None:
					found_directional_word = 0 # set this to 0 so it can search other Gleason scores in the report
					if primary_grade + secondary_grade <= 6:
						overall_grade_group = 1
					elif primary_grade + secondary_grade == 7:
						if primary_grade == 3:
							overall_grade_group = 2
						elif primary_grade == 4:
							overall_grade_group = 3
						else: # 7/10 with no primary and seoncdary info
							overall_grade_group = 2
					elif primary_grade + secondary_grade == 8:
						overall_grade_group = 4
					elif primary_grade + secondary_grade > 8:
						overall_grade_group = 5
					else:
						overall_grade_group = None

					# update the grade group lists
					if overall_grade_group is not None:
						primary_grade_list_local.append(primary_grade)
						secondary_grade_list_local.append(secondary_grade)
						overall_grade_list_local.append(overall_grade_group)
						overall_gs_list_local.append(primary_grade + secondary_grade)
				else:
					overall_grade_group = None

				if num_pos_core is not None:
					num_pos_core_list_local.append(num_pos_core)
					num_total_core_list_local.append(num_total_core)

				if len(max_core_involve) > 0:# is not None:
					max_core_involve_list_local.append(max_core_involve)
				
			# primary, secondary, overall grade groups
			if len(primary_grade_list_local) > 0:
				# print('Report Number : ', report_number)	
				# print('text_oi : ')
				# print(text_oi)
				# print('overall_grade_list_local : ', overall_grade_list_local)
				# breakpoint()
				primary_grade_list.append(primary_grade_list_local)
				secondary_grade_list.append(secondary_grade_list_local)
				overall_grade_list.append(overall_grade_list_local)	
				if len(set(overall_grade_list_local)) > 0: # just get average of local overall grades
					# overall_grade_merged_list.append(overall_grade_list_local[0])
					overall_grade_merged_list.append(np.max(overall_grade_list_local))
				else:
					overall_grade_merged_list.append(None)
				overall_gs_list.append(overall_gs_list_local)
			else:
				primary_grade_list.append(None)
				secondary_grade_list.append(None)
				overall_grade_list.append(None)	
				overall_gs_list.append(None)
				overall_grade_merged_list.append(None)

			# num cores
			if len(num_pos_core_list_local) > 0: # if no postiive nodes identified
				num_pos_cores_list.append(num_pos_core_list_local)
				num_pos_cores_sum_list.append(np.sum(num_pos_core_list_local))
				num_total_core_list.append(num_total_core_list_local)
				num_total_core_sum_list.append(np.sum(num_total_core_list_local))
				# if np.sum(num_total_core_list_local) == 30:
				# 	print('num_pos_core')
				# 	print(num_total_core_list_local)
				# 	print(text_oi)
				# 	breakpoint()
				# if np.sum(num_pos_core_list_local) == 17:
				# 	print('num_pos_core')
				# 	print(num_pos_core_list_local)
				# 	print(text_oi)
				# 	print(raw_gleason_cores_context)
				# 	breakpoint()
				# if np.sum(num_pos_core_list_local) == 18:
				# 	print('num_pos_core')
				# 	print(num_pos_core_list_local)
				# 	print(text_oi)
				# 	print(raw_gleason_cores_context)
				# 	breakpoint()
				# if np.sum(num_pos_core_list_local) == 16:
				# 	print('num_pos_core')
				# 	print(num_pos_core_list_local)
				# 	print(text_oi)
				# 	print(raw_gleason_cores_context)
				# 	breakpoint()
				# if np.sum(num_pos_core_list_local) == 15:
				# 	print('num_pos_core')
				# 	print(num_pos_core_list_local)
				# 	print(text_oi)
				# 	print(raw_gleason_cores_context)
				# 	breakpoint()
				# if np.sum(num_total_core_list_local) == 20:# == num_total_core:
				# 	print('num_total_core')
				# 	print(num_total_core_list_local)
				# 	print(text_oi)
				# 	breakpoint()
				# if np.sum(num_total_core_list_local) == 40:# == num_total_core:
				# 	print('num_total_core')
				# 	print(num_total_core_list_local)
				# 	print(text_oi)
				# 	breakpoint()
				# if np.sum(num_total_core_list_local) == 80:# == num_total_core:
				# 	print('num_total_core')
				# 	print(num_total_core_list_local)
				# 	print(text_oi)
				# 	breakpoint()
			else:
				num_pos_cores_list.append(None)
				num_total_core_list.append(None)
				num_pos_cores_sum_list.append(None)
				num_total_core_sum_list.append(None)

			# max core involvement
			if len(max_core_involve_list_local) > 0:
				max_core_involve_list.append(max_core_involve_list_local)
			else:
				max_core_involve_list.append(None)

			# histology
			if adenocarcinoma_ind:
				adenocarcinoma.append(1)
			else:
				adenocarcinoma.append(0)

			if small_cell_carc_ind and not adenocarcinoma_ind:
				small_cell_carc.append(1)
			else:
				small_cell_carc.append(0)
				
			if neuroendocrine_carc_ind and not adenocarcinoma_ind:
				neuroendocrine_carc.append(1)
			else:
				neuroendocrine_carc.append(0)
				
			# tumor stage
			if not tumor_stage_ind:
				tumor_stage_list.append(None)
			# margin
			if not margin_ind:
				margin_list.append(None)

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_pos_cores'] = num_pos_cores_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'primary_grade'] = primary_grade_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'secondary_grade'] = secondary_grade_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'overall_grade_group'] = overall_grade_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'overall_gs'] = overall_gs_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'overall_grade_merged'] = overall_grade_merged_list

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'benign'] = benign_list

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_pos_cores'] = num_pos_cores_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_total_core'] = num_total_core_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_pos_cores_sum'] = num_pos_cores_sum_list
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'num_total_core_sum'] = num_total_core_sum_list

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'max_core_involve'] = max_core_involve_list
		# manual gleason score fix (primary, secondary) : # (3,7) => (3,4) Report MS04R67776, (3,6) => (3,3) Report S0009791E
		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number == 'MS04R67776', 'secondary_grade'] = 4
		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number == 'MS04R67776', 'overall_grade_group'] = 2

		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number == 'S0009791E', 'secondary_grade'] = 3
		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number == 'S0009791E', 'overall_grade_group'] = 1

		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'radical_prostatectomy'] = radical_pros

		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'small_cell_carc'] = small_cell_carc
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'neuroendocrine_carc'] = neuroendocrine_carc
		df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'adenocarcinoma'] = adenocarcinoma
		# breakpoint()

		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'margin'] = margin_list
		# # manually fix margin based on Madhur's audit
		# positive_margin_reports = ['S94W17361', 'S95R18902', 'MS09K38970', 'BS09R32983', 'S16-17016']
		# negative_margin_reports = ['S9830284L', 'S0133510G', 'S0020803D', 'BS07T08214', 'MS09G67481', 'MS08K66009', 'S9818403C']
		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number.isin(positive_margin_reports), 'margin'] = 1
		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number.isin(negative_margin_reports), 'margin'] = 0

		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'pT_stage'] = [None if stage is None else 'p' + stage for stage in tumor_stage_list]
		# pt_stage_list = []
		# remove_chars = ['o', 'x', 's']
		# for pt_stage in df_pathology_with_biopsy.pT_stage.values:
		# 	if pt_stage is None:
		# 		pt_stage_list.append(None)
		# 		continue
		# 	if pt_stage == 'pt30':
		# 		pt_stage_list.append('pt3')
		# 		continue
		# 	elif pt_stage == 'pt20':
		# 		pt_stage_list.append('pt2')
		# 		continue
		# 	elif pt_stage == 'pt31': # manual override in accesion number : bs13j21799
		# 		pt_stage_list.append('pt3a')
		# 		continue
		# 	elif pt_stage == 'pt25': 
		# 		pt_stage_list.append('pt2')
		# 		continue
		# 	elif pt_stage == 'pt26': 
		# 		pt_stage_list.append('pt2')
		# 		continue
		# 	pt_stage_filtered = pt_stage
		# 	for char in remove_chars:
		# 		pt_stage_filtered = pt_stage_filtered.replace(char, '')
		# 	pt_stage_list.append(pt_stage_filtered)
		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'pT_stage'] = pt_stage_list
		# df_pathology_with_biopsy = get_tumor_stage_info_via_context(df_pathology_with_biopsy)
		# # test and statistics
		# df_pathology_with_biopsy_comparison = df_pathology_with_biopsy.dropna(subset = ['pT_stage_context'])
		# df_pathology_with_biopsy_comparison = df_pathology_with_biopsy_comparison.dropna(subset = ['pT_stage'])
		# df_compare = df_pathology_with_biopsy_comparison[['Report_Number', 'pT_stage', 'pT_stage_context']]
		# print('\n')
		# print('Discrepancies between the numeric and context-based methods among the reports which have the both information available')
		# count = 0
		# for entry in df_compare.values:
		# 	compare_len = len(entry[2])
		# 	if entry[2] != entry[1][0:compare_len]:
		# 		print(entry)
		# 		count = count + 1
		# print('AJCC numeric-based pT stage vs. context-based pT stage')
		# print('accuracy : ', (len(df_compare) - count)/len(df_compare))
		# print('\n')
		# print('Number of reports in each pT stage : ')
		# # combining numeric based and context based pT stage 
		# pt_stage_combined = []; context_based_counter = 0; numeric_based_counter = 0;
		# for pt_stage_numeric, pT_stage_context in zip(df_pathology_with_biopsy.pT_stage.values, df_pathology_with_biopsy.pT_stage_context.values):
		# 	if pt_stage_numeric is not None:
		# 		pt_stage_combined.append(pt_stage_numeric)
		# 		numeric_based_counter += 1
		# 	elif pT_stage_context is not None:
		# 		pt_stage_combined.append(pT_stage_context)
		# 		context_based_counter += 1
		# 	else:
		# 		pt_stage_combined.append(None)
		# df_pathology_with_biopsy.loc[df_pathology_with_biopsy.index, 'pT_stage_combined'] = pt_stage_combined
		# df_count_stats = pd.value_counts(df_pathology_with_biopsy.pT_stage_combined.values, sort = True)
		# print(df_count_stats)
		# print('# of numeric-based pT stages : ', numeric_based_counter)
		# print('# of context-based pT stages : ', context_based_counter)
		# print('Among both methods, {0} is the number of reports with hard-coded pT stage sentences'.format(pt_stage_counter_hard_coded))
		# print('total = ', sum(df_count_stats.values))
		# print('\n')

	if not load_prcossed:
		primary_grade_all_list = []
		for val in df_pathology_with_biopsy.primary_grade.values:
			if val is not None:
				primary_grade_all_list += val
		df_count_stats_primary_grade = pd.value_counts(primary_grade_all_list, sort = True)		
		print('Counts per each primary grade: ')
		print(df_count_stats_primary_grade)
		print('\n')

		secondary_grade_all_list = []
		for val in df_pathology_with_biopsy.secondary_grade.values:
			if val is not None:
				secondary_grade_all_list += val
		df_count_stats_secondary_grade = pd.value_counts(secondary_grade_all_list, sort = True)		
		print('Counts per each secondary grade: ')
		print(df_count_stats_secondary_grade)
		print('\n')

		# outlier filtering
		max_core_involve_all_list = []; outlier_report_list = []
		mci_outliers = [890, 0, 800, 320]
		for val, report_number_mci in zip(df_pathology_with_biopsy.max_core_involve.values, df_pathology_with_biopsy.Report_Number.values):
			if val is not None:
				for val_sub in val:
					if any([val_outlier in val_sub for val_outlier in mci_outliers]):
						outlier_report_list.append(report_number_mci)
					else:
						max_core_involve_all_list += val_sub
		df_count_stats_max_core_involve = pd.value_counts(max_core_involve_all_list, sort = True)		
		print('Counts per each max core involvement : ')
		print(df_count_stats_max_core_involve)
		print('\n')
		# remove outliers
		reports_mci_oi = set(df_pathology_with_biopsy.Report_Number.values) - set(outlier_report_list)
		# breakpoint()
		df_pathology_with_biopsy = df_pathology_with_biopsy.loc[df_pathology_with_biopsy.Report_Number.isin(reports_mci_oi)]

		# benign stats
		df_count_stats_benign = pd.value_counts(df_pathology_with_biopsy.benign.values, sort = True)
		print('Number of reports with benign indicator 1 (postiive) or 0 (negative) : ')
		print(df_count_stats_benign)
		# print('Among those, {0} is the number of reports with hard-coded margin indicator sentences'.format(margin_counter_hard_coded))
		print('total = ', sum(df_count_stats_benign.values))
		print('\n')
		
		# overall grade stats
		df_count_stats_overall_grade = pd.value_counts(df_pathology_with_biopsy.overall_grade_merged.values, sort = True)

		df_count_stats_num_pos_cores = pd.value_counts(df_pathology_with_biopsy.num_pos_cores_sum.values, sort = True)
		# outlier filtering
		num_pos_sum_outliers = [20, 30]
		num_pos_cores_oi = set(df_count_stats_num_pos_cores.loc[df_count_stats_num_pos_cores.values >= 9].index) | set([None]) - set(num_pos_sum_outliers)
		df_pathology_with_biopsy = df_pathology_with_biopsy.loc[df_pathology_with_biopsy.num_pos_cores_sum.isin(num_pos_cores_oi)]

		df_count_stats_num_total_cores = pd.value_counts(df_pathology_with_biopsy.num_total_core_sum.values, sort = True)
		# outlier filtering
		num_total_core_sum_outliers = [30]
		num_total_cores_oi = set(df_count_stats_num_total_cores.loc[df_count_stats_num_total_cores.values >= 9].index) | set([None]) - set(num_total_core_sum_outliers)
		df_pathology_with_biopsy = df_pathology_with_biopsy.loc[df_pathology_with_biopsy.num_total_core_sum.isin(num_total_cores_oi)]

		# get new stats
		df_count_stats_num_pos_cores = pd.value_counts(df_pathology_with_biopsy.num_pos_cores_sum.values, sort = True)
		df_count_stats_num_total_cores = pd.value_counts(df_pathology_with_biopsy.num_total_core_sum.values, sort = True)
		print('Number of reports per each # of positive cores : ')
		print(df_count_stats_num_pos_cores)
		print('Number of reports per each # of total cores : ')
		print(df_count_stats_num_total_cores)
		print('Number of reports per each overall grade : ')
		print(df_count_stats_overall_grade)
		print('total = ', df_count_stats_overall_grade.sum())
		# breakpoint()
		# df_pathology_with_biopsy_oi = df_pathology_with_biopsy.copy()#.loc[df_pathology_with_biopsy.radical_prostatectomy == 1]
		print('Number of unique biopsy reports : ', len(df_pathology_with_biopsy))

		print('Number of unique biopsy reports per center : ')
		print(pd.value_counts(df_pathology_with_biopsy.MRN_Type.values, sort = True))

		print('Number of unique patients : ', len(set(df_pathology_with_biopsy.EMPI.values)))

		df_pathology_with_biopsy_oi_overall_grade = df_pathology_with_biopsy.dropna(subset = ['overall_grade_group'])
		print('Number of unique patients with overall grade : ', len(set(df_pathology_with_biopsy_oi_overall_grade.EMPI.values)))

		df_pathology_with_biopsy_oi_with_num_cores_involved = df_pathology_with_biopsy_oi_overall_grade.dropna(subset = ['num_pos_cores'])
		df_pathology_with_biopsy_oi_with_num_cores_involved_wo_grade = df_pathology_with_biopsy.dropna(subset = ['num_pos_cores'])
		print('Number of unique patients with overall grade, and num. of positive and total cores : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved.EMPI.values)))
		print('Number of unique patients with num. of positive and total cores : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_wo_grade.EMPI.values)))


		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve = df_pathology_with_biopsy_oi_with_num_cores_involved.dropna(subset = ['max_core_involve'])
		print('Number of unique patients with overall grade, num. of positive and total cores, and maximum core involvement : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve.EMPI.values)))
		
		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign = df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve.dropna(subset = ['benign'])
		print('Number of unique patients with overall grade, num. of positive and total cores, maximum core involvement, and benign : ', len(set(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.EMPI.values)))
		print('\n')
		print('\n')
		# breakpoint()
	continue_ = True; sample_rp = True; specific_report = False
	if load_prcossed:
		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign = df_pathology_with_biopsy
	
	df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.drop(columns = ['Unnamed: 0'], inplace = True)
	df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.index = np.arange(len(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign))
	df_biopsy_wo_duplicate_processed = featurize_biopsy_df(df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign)
	df_biopsy_wo_duplicate_processed.to_csv('../data/df_biopsy_wo_duplicate_processed.csv')
	breakpoint()
	while continue_:
		# df_pathology_with_gleason
		print('\n')
		print('Extracted info : ')
		# if sample_rp:
		# 	sampled_entry = df_pathology_with_biopsy_oi_with_pt_stage_margin.loc[df_pathology_with_biopsy_oi_with_pt_stage_margin.radical_prostatectomy == 1].sample(n = 1)
		# else:
		if specific_report:
			sampled_entry = df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.loc[df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.Report_Number == report_number_to_query]
		else:
			sampled_entry = df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.sample(n = 1)
		sampled_entry_extracted_oi = sampled_entry[['Report_Number', 'overall_gs', 'overall_grade_group', 'overall_grade_merged', 'num_pos_cores', 'num_pos_cores_sum', 'num_total_core', 'num_total_core_sum', 'max_core_involve', 'benign', 'small_cell_carc', 'neuroendocrine_carc', 'adenocarcinoma']]
		print(sampled_entry_extracted_oi.iloc[0])
		print('\n')
		print('Pathology report : ')
		print(sampled_entry.Report_Text.values[0])
		print('\n')
		continue_ = input('Would you like to continue? (Y/N) : ') == 'Y'
		specific_report = input('Want to look at specific report? (Y/N) : ') == 'Y'
		if specific_report:
			report_number_to_query = input('Type report number : ')
			if report_number_to_query not in df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.Report_Number.values:
				print('Report not found!')
				print('Randomly sampling a report...')
				specific_report = False
				continue
	
	export = input('Export the pathology report datafame with the extracted info in the current directory? (Y/N) : ')
	if export not in {'Y', 'N'}:
		raise KeyError('Input Y or N only')
	if export == 'Y':
		df_pathology_with_biopsy.to_csv('df_pathology_with_biopsy_params_extracted.tsv', sep = '\t', index = True)
		df_pathology_with_biopsy_oi_with_num_cores_involved_max_core_involve_benign.to_csv('df_pathology_with_biopsy_params_extracted_all_feats.tsv', sep = '\t', index = True)
		print('\n')
		print('The processed file exported to ' + os.getcwd() + '/df_pathology_with_extracted_info.tsv')
		print('\n')
	return

if __name__ == "__main__":
	explain = """
	Written by Intae Moon
	Date : May 3rd, 2021

	Description :
	The following script extract 
	1. Overall grade group
	2. Benign indicator
	3. Number of positive cores
	4. Total number of cores
	5. Maximum core involvement

	See : https://docs.google.com/document/d/1pRA2XcAxjbqbji8WHFfJD2bvwgo-EgryGLNRiq-gYT4/edit?usp=sharing for more details
	"""
	print("\n")
	print(explain)
	print("\n")
	input('Type any button if you would like to continue')
	main()
