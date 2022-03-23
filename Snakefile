configfile: "config.yaml"
from scripts.utils import *            
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]
Sample_dict = Get_Sample_dict(data_dir)
Projects = Sample_dict.keys()
Samples = [sample for project  in Projects for sample in Sample_dict[project].values()]
paths = Get_file_paths(Sample_dict, data_dir)

#-------------------------------------------------------------------------------------------------------------------
# 0.3 specify target rules
rule all:
    input:
        expand("../data/copynumber/ACE/{sample}/2N/fits.txt", sample = Samples)
        

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
        errorgraph =  "../data/copynumber/ACE/{sample}/2N/errorgraph.svg"
    params:
        ploidies = config['ACE']['ploidies'],
        method = config['ACE']["method"],
        penalty = config['ACE']['penalty']
    conda:
        "envs/ACE.yaml"
    script:
        "scripts/ACE.R"
    
        
