#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# utils.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Utils code for TCGA_purity snakemake
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1)
#
# History:
#  09-03-2022: Initialize script, Add: Get_Sample_dict, Get_file_paths, Migrate_files
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
import pandas as pd
import glob
from shutil import copyfile
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1.1 Specify functions
#-------------------------------------------------------------------------------
def Get_Sample_dict(data_dir):
    """ 
    Obtain dictionary with file information per sample
    :param str data_dir: data directory where sampletables are stored
    :returns: dict Sample_dict: nested dict per project with sample IDs
    """
    Sample_dict = dict()
    # read CNA_sampletable
    CNA_sampletable = pd.read_csv(data_dir + "sampletable/CNA_Sampletable.tsv", sep = "\t")
    # filter out Normals
    CNA_sampletable = CNA_sampletable[CNA_sampletable['sample_type'] == "Primary Tumor"] 
    # Obtain TCGA sample id: split on 3rd "-"
    CNA_sampletable['TCGA_sample'] = CNA_sampletable['cases'].apply(lambda x: "-".join(x.split("-")[:3]))
    # iterate over rows
    for project,hash_sample,TCGA_sample in zip(CNA_sampletable.project,CNA_sampletable.id,CNA_sampletable.TCGA_sample):
        # create nested directory
        if project not in Sample_dict:
            Sample_dict[project] = dict()
        # add samples to dict
        Sample_dict[project][hash_sample] = TCGA_sample
    return Sample_dict

def Get_file_paths(Sample_dict, data_dir):
    """ 
    Obtain list with all input files
    :param dict Sample_dict:  dictionary with file information per sample
    :param str data_dir: data directory where sampletables are stored
    :returns: list paths: list of paths to input files
    """
    paths = []
    # iterate over projects
    for project,id_dict in Sample_dict.items():
        # obtain subtype
        subtype = project.split("-")[-1]
        # iterate over file ids
        for id in id_dict:
            # glob the filename
            path = glob.glob(f"{data_dir}{subtype}/{project}/harmonized/Copy_Number_Variation/Masked_Copy_Number_Segment/{id}/*.txt")[0]
            paths.append(path)
    return paths
            
        
def Migrate_files(Sample_dict,paths, data_dir):
    """ 
    :param dict Sample_dict:  dictionary with file information per sample
    :param list paths: list of paths to input files
    :param str data_dir: data directory where sampletables are stored
    :returns: list paths: list of paths to input files
    """
    for path in paths:
        project = path.split("/")[3]
        id = path.split("/")[-2]
        TCGA_sample = Sample_dict[project][id]
        copyfile(path,f"{data_dir}/copynumber/segments/{TCGA_sample}_segments.txt")
        
        
        
