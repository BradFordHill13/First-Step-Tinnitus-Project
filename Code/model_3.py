#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 17:12:41 2020

@author: matthieusimeoni
"""
import pandas as pd
import numpy as np
import sys
            
# read samples file: create a dataframe that contains uniquely the participant ids
def read_samples(sample_file):
    samples = pd.read_csv(sample_file, sep=" ", dtype={'ID_1': int,'ID_2': int,'missing': int,'sex': str})
    # drop test value at index 0: "0 0 0 D"
    samples.drop(samples.index[0], inplace=True) 
    # keep non-missing samples only
    samples = samples[samples['missing'] == 0]
    # drop column that are not needed (missing and sex)
    samples.drop(samples.columns[[1,2,3]], axis=1, inplace=True)
    # rename the one column that is needed: eid
    samples.rename(columns = {'ID_1':'eid'}, inplace = True)
    samples.set_index("eid", inplace = True)
    return samples


def convert_tinnitus_legend(raw_value): # do you have TINNITUS?
    if raw_value==0: return 0 # no
    elif raw_value==11: return 4 # yes, always
    elif raw_value==12: return 3 # yes, often
    elif raw_value==13: return 2 # yes, sometimes
    elif raw_value==14: return 0 # not now, but I had it in the past = no
    else: return np.nan


def convert_tinnitus_legend_score(raw_value): # do you have TINNITUS?
    if raw_value==4: return -0.5 # no
    elif raw_value==11: return 1.0 # yes, always
    elif raw_value==12: return 0.5 # yes, often
    elif raw_value==13: return 0.0 # yes, sometimes
    else: return 0.0


def assign_tinnitus_phenotype(row):
    # retrieve tinnitus score (pheno id 4803)
    n = 0 # number of values
    selected_values_1 = []
    if not np.isnan(row["4803-0.0"]): 
        tinnitus_status_0 = convert_tinnitus_legend(row["4803-0.0"]) 
        n += 1
        selected_values_1.append(tinnitus_status_0)
    if not np.isnan(row["4803-1.0"]): 
        tinnitus_status_1 = convert_tinnitus_legend(row["4803-1.0"]) 
        n += 1
        selected_values_1.append(tinnitus_status_1)
    if not np.isnan(row["4803-2.0"]): 
        tinnitus_status_2 = convert_tinnitus_legend(row["4803-2.0"]) 
        n += 1
        selected_values_1.append(tinnitus_status_2)
    if not np.isnan(row["4803-3.0"]): 
        tinnitus_status_3 = convert_tinnitus_legend(row["4803-3.0"]) 
        n += 1
        selected_values_1.append(tinnitus_status_3)
    
    for index, i in enumerate(selected_values_1): # NaN management 
        if np.isnan(i): # if there is NaN, we remove it
            list_no_nan = selected_values_1
            list_no_nan.remove(selected_values_1[index])
            n -= 1 # because we counted 1 before for NaN
        else: # else if there is no NaN
            list_no_nan = selected_values_1
    
    selected_values_1 = list_no_nan
    
    def all_same(items): # to control in case there is [0, 0, 0, 0]
        if len(items) > 1:
            return all(x == items[0] for x in items)
    
    list_no_zero = []
    for index, i in enumerate(selected_values_1): # 0 management
        if all_same(selected_values_1) == True and i == 0:
            list_no_zero = selected_values_1
            break
        elif i == 0 and len(selected_values_1) > 1: # if there is 0 and others numbers, we remove 0
            list_no_zero = selected_values_1
            list_no_zero.remove(selected_values_1[index]) # we remove 0 because 0 do not average
            n -=1 # because we counted 1 before for 0
        else:
            list_no_zero = selected_values_1
    
    selected_values_1 = list_no_zero
    
    if n != 0:
        mean_selected_values = sum(selected_values_1)/n # if we add round() then 2.5 will be rounded up to 2
    else:
        mean_selected_values = np.nan
    
    # assignment of severity level
    selected_values_2 = []
    if not np.isnan(row["4814-0.0"]):
        tinnitus_severity_score_0 = convert_tinnitus_legend_score(row["4814-0.0"])
        selected_values_2.append(tinnitus_severity_score_0)
    if not np.isnan(row["4814-1.0"]):
        tinnitus_severity_score_1 = convert_tinnitus_legend_score(row["4814-1.0"])
        selected_values_2.append(tinnitus_severity_score_1)
    if not np.isnan(row["4814-2.0"]):
        tinnitus_severity_score_2 = convert_tinnitus_legend_score(row["4814-2.0"])
        selected_values_2.append(tinnitus_severity_score_2)
    if not np.isnan(row["4814-3.0"]):
        tinnitus_severity_score_3 = convert_tinnitus_legend_score(row["4814-3.0"])
        selected_values_2.append(tinnitus_severity_score_3)
    if len(selected_values_2) == 0: # if it has only missing values, then the score is 0
        selected_values_2.append(0) # add to the list and do not transform to int because by using sum() on int, the latter is not iterable
    
    if len(selected_values_2) == 1 and selected_values_2[0] == 0: # if it's 0 
        tinnitus_severity_score = 0
    else:
        tinnitus_severity_score = sum(selected_values_2) # if we add round() then 2.5 will be rounded up to 2

    # return the final score
    tinnitus_status = mean_selected_values + tinnitus_severity_score
    
    if tinnitus_status < 0: # if the status is less than 0 then it's a control
        tinnitus_status = 0
    return tinnitus_status


def tinnitus_extraction_to_phenofile(output_file, sample_file, extraction_file):
    
    # create an empty dataframe to store phenofile output
    # as index, use the eid (i.e. the UKBB participant ids) from sample file (to respect UKBB ordering)
    output = read_samples(sample_file)
    output['tinnitus'] = np.nan # add tinnitus phenotype column, to be populated
    
    # read data from UKBB tinnitus extraction file 
    # (4803-x.y: tinnitus status, 4814-x.y: severity level)
    tinnitus_dtype={'eid': int,\
                    '4803-0.0': float,'4803-1.0': float,'4803-2.0': float, '4803-3.0': float,\
                    '4814-0.0': float,'4814-1.0': float,'4814-2.0': float,'4814-3.0': float}
    tinnitus_extraction = pd.read_csv(extraction_file, sep=",", dtype=tinnitus_dtype)
    tinnitus_extraction.set_index('eid', inplace = True) # use eid as index

    # loop over extracted eids (i.e. UKBB participant ids)
    not_in_sample=[]
    counter=0
    print("Please wait while UK Biobank extraction is processed...")
    for index, _ in output.iterrows():
        eid = index
        try: # eid might not be present in tinnitus_phenotype
            tinnitus_extraction_row = tinnitus_extraction.loc[eid]
            output.loc[eid] = assign_tinnitus_phenotype(tinnitus_extraction_row)
        except KeyError:
            eid_not_in_sample = str(eid)
            not_in_sample.append(eid_not_in_sample)
        counter+=1 # keep track of number of processed rows
        if counter%10000==0: sys.stdout.write('.') # print out once in a while 
     
    valid = output.iloc[:,0].notna().sum() # count number of non-NA results
    n_cases = len(output[(output["tinnitus"]>0)])
    n_controls = len(output[(output["tinnitus"]==0)])
    print("\nPhenotype file has been built:")
    print(valid,"entries were included:",n_cases,"cases and ",n_controls,"controls")
    print(len(not_in_sample), "entries were omitted (not in sample file)")
    output.fillna(-999, inplace=True) # replace NAs according to BGENIE convention
    output.replace(-1, -999, inplace=True) # -1 is also equivalent to an NA, for tinnitus
    # output
    output.to_csv(output_file, sep=' ', index = False)
    omitted_file = output_file[:-4] + ".NotInSamplefile.csv"
    with open(omitted_file, "w") as outfile:
        outfile.write("\n".join(not_in_sample))


def main():
    ###########################################################################
    # CONFIG SEVERITY
    # "any"
    # or a string containing a combination of severe, moderate, slight, e.g.: "severe-moderate"
    ###########################################################################
    output_file = "./output/phenofile_tinnitus_model_3.csv" # "./test_output/phenofile_tinnitus_model_3_test.csv"
    sample_file = "input/ukb43805_imp_chr1_v3_s487297.sample" # list of UKBB participants # "test_input/samples_order.sample"
    extraction_file = "input/ukb42625.csv" # "test_input/model_1_input.csv"
    print("Combining tinnitus data (from UK Biobank extraction) to generate phenotype file for GWAS")
    print("Using sample file (list of participant):", sample_file)
    print("Severity:\n")
    tinnitus_extraction_to_phenofile(output_file, sample_file, extraction_file)
    print("DONE")
  
if __name__== "__main__":
    main()