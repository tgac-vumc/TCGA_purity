configfile: "config.yaml"
from scripts.utils import *            
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]
Sample_dict = Get_Sample_dict(data_dir)
Projects = Sample_dict.keys()
Samples = set([sample for project  in Projects for sample in Sample_dict[project].values()])
paths = Get_file_paths(Sample_dict, data_dir)
Purity_Measures = ['ESTIMATE','ABSOLUTE','LUMP','IHC', 'ACE']
#-------------------------------------------------------------------------------------------------------------------
# 0.3 specify target rules
rule all:
    input:
        "../plots/purity/Scatterplots_purity.svg"
        

#++++++++++++++++++++++++++++++++++++++++++++++++ 1 PREPARE DATA  +++++++++++++++++++++++++++++++++++++++++++++++++
# Migrate segmentfiles
rule Migrate_files:
    input:
        paths
    output:
        dynamic("../data/copynumber/segments/{sample}_segments.txt")
    run:
        Migrate_files(Sample_dict, input, data_dir)


#++++++++++++++++++++++++++++++++++++++++++ 2 ESTIMATE PURITY WITH ACE  +++++++++++++++++++++++++++++++++++++++++++
# Run ACE using a custom ACE script

rule ACE:
    input:
        segments = "../data/copynumber/segments/{sample}_segments.txt"
    output:
        fit = "../data/copynumber/ACE/{sample}/2N/fits.txt",
        errorgraph = "../data/copynumber/ACE/{sample}/2N/errorgraph.svg"
    params:
        binsize = config['ACE']['binsize'],
        ploidies = config['ACE']['ploidies'],
        method = config['ACE']["method"],
        penalty = config['ACE']['penalty']
    conda:
        "envs/ACE.yaml"
    script:
        "scripts/ACE.R"
    
        
#++++++++++++++++++++++++++++++++++++++++++ 3 COMPARE PURITY MEASURES   +++++++++++++++++++++++++++++++++++++++++++
# Plot scatterplots of purity estimates of different measures

rule Compare_purity_measures:
    input:
        ACE_fits = expand("../data/copynumber/ACE/{sample}/2N/fits.txt", sample = Samples),
        Purities = "../data/purity/Tumor_purities.tsv"
    output:
        Scatterplots = "../plots/purity/Scatterplots_purity.svg",
    conda:
        "envs/R.yaml"
    script:
        "scripts/Compare_purity_measures.R"


#+++++++++++++++++++++++++++++++++++++ 4 ESTIMATE CONSENSUS PURITY WITH ACE  ++++++++++++++++++++++++++++++++++++++
# Estimate consensus purity including ACE (instead of ABSOLUTE)
"""
rule Consensus_purity_estimation:
    input:
        segments = "../data/copynumber/segments/{sample}_segments.txt"
    output:
        fit = "../data/copynumber/ACE/{sample}/2N/fits.txt",
        errorgraph =  "../data/copynumber/ACE/{sample}/2N/errorgraph.svg"
    params:
        ploidies = config['ACE']['ploidies'],
        method = config['ACE']["method"],
        penalty = config['ACE']['penalty']
    conda:
        "envs/ACE.yaml"
    script:
        "scripts/ACE.R"
"""
