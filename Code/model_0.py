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
    elif raw_value==11: return 3 # yes, always
    elif raw_value==12: return 2 # yes, often
    elif raw_value==13: return 1 # yes, sometimes
    elif raw_value==14: return np.nan # not now, but I had it in the past
    else: return np.nan

def assign_tinnitus_phenotype(row, severity_level):
    # retrieve tinnitus score (pheno id 4803)
    tinnitus_status = convert_tinnitus_legend(row["4803-0.0"]) # if first assessment is valid, keep it
    if np.isnan(tinnitus_status):
        tinnitus_status = convert_tinnitus_legend(row["4803-1.0"]) # if not, keep going...
    if np.isnan(tinnitus_status):
        tinnitus_status = convert_tinnitus_legend(row["4803-2.0"])
    if np.isnan(tinnitus_status):
        tinnitus_status = convert_tinnitus_legend(row["4803-3.0"])
    
    # return score irrespective of severity, if level set to "any"...
    if severity_level=="any": return tinnitus_status
    # return score irrespective of severity, if subject did not experience tinnitus
    if tinnitus_status==0: return tinnitus_status
    
    # ...otherwise, retrieve tinnitus severity score (pheno id 4803)
    tinnitus_severity_score = row["4814-0.0"] # if first assessment is valid, keep it
    if np.isnan(tinnitus_severity_score):
        tinnitus_severity_score = row["4814-1.0"] # if not, keep going...
    if np.isnan(tinnitus_severity_score):
        tinnitus_severity_score = row["4814-2.0"]
    if np.isnan(tinnitus_severity_score):
        tinnitus_severity_score = row["4814-3.0"]
    
    # if retrieved severity score matches user-defined severity level, return tinnitus_status
    if "severe" in severity_level and tinnitus_severity_score==11:
        return tinnitus_status
    elif "moderate" in severity_level and tinnitus_severity_score==12:
        return tinnitus_status
    elif "slight" in severity_level and tinnitus_severity_score==13:
        return tinnitus_status
    else: return np.nan # ...otherwise, return nan as tinnitus phenotype


def tinnitus_extraction_to_phenofile(output_file, sample_file, extraction_file, severity_level):
    
    # create an empty dataframe to store phenofile output
    # as index, use the eid (i.e. the UKBB participant ids) from sample file (to respect UKBB ordering)
    output = read_samples(sample_file)
    output['FID'] = np.nan # add FID according to REGENIE convention
    output['IID'] = np.nan # add FID according to REGENIE convention
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
            output.loc[eid] = assign_tinnitus_phenotype(tinnitus_extraction_row, severity_level)
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
    output['tinnitus'] = output['tinnitus'].astype(int) # float to int for all the values
    output['tinnitus'] = output.replace(-999, np.nan, inplace=False) # replace in NAs according to REGENIE convention
    # output
    output['FID'] = output.index
    output['IID'] = output.index
    # Add interger and NaN in tinnitus column using str (int doesn't work with NaN because NaN is float)
    output['tinnitus'] = output['tinnitus'].fillna('NA')
    output['tinnitus'] = output['tinnitus'].astype(str)
    output['tinnitus'] = output['tinnitus'].str.split('.')
    output['tinnitus'] = output['tinnitus'].str[0]
    output.to_csv(output_file, sep=' ', index = False)
    omitted_file = output_file[:-4] + ".NotInSamplefile.csv"
    with open(omitted_file, "w") as outfile:
        outfile.write("\n".join(not_in_sample))


def main():
    ###########################################################################
    # CONFIG SEVERITY
    # "any"
    # or a string containing a combination of severe, moderate, slight, e.g.: "severe-moderate"
    severity_level = "any"
    ###########################################################################
    output_file = "./output/phenofile_tinnitus_LOG_model_0_1and2and3_"+severity_level+".csv"
    #output_file ="./test_output/phenofile_tinnitus_LOG_model_0_1and2and3_"+severity_level+".csv"
    sample_file = "input/ukb43805_imp_chr1_v3_s487297.sample" # list of UKBB participants
    #sample_file = "test_input/samples_order.sample"
    extraction_file = "input/ukb42625.csv"
    #extraction_file = "test_input/model_0_input.csv"
    print("Combining tinnitus data (from UK Biobank extraction) to generate phenotype file for GWAS")
    print("Using sample file (list of participant):", sample_file)
    print("Severity:",severity_level,"\n")
    tinnitus_extraction_to_phenofile(output_file, sample_file, extraction_file, severity_level)
    print("DONE")
  
if __name__== "__main__":
    main()
