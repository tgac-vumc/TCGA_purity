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
        expand("../plots/purity/Scatterplots_purity_{ploidy}.svg", ploidy = config['ACE']['ploidies'])

#++++++++++++++++++++++++++++++++++++++++++++++++ 1 PREPARE DATA  +++++++++++++++++++++++++++++++++++++++++++++++++
# Migrate segmentfile
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
        fit = "../data/copynumber/ACE/{sample}/{ploidy}/fits.txt",
        errorgraph = "../data/copynumber/ACE/{sample}/{ploidy}/errorgraph.svg"
    params:
        method = config['ACE']["method"],
        penalty = config['ACE']['penalty'],
        penploidy = config['ACE']['penploidy']
    conda:
        "envs/ACE.yaml"
    script:
        "scripts/ACE.R"
    
#++++++++++++++++++++++++++++++++++++++++++ 3 COMPARE PURITY MEASURES   +++++++++++++++++++++++++++++++++++++++++++
# Plot scatterplots of purity estimates of different measures

rule Compare_purity_measures:
    input:
        ACE_fits = lambda wildcards: expand("../data/copynumber/ACE/{sample}/"+ wildcards.ploidy + "/fits.txt", sample = Samples),
        Purities = "../data/purity/Tumor_purities.tsv"
    output:
        Scatterplots = "../plots/purity/Scatterplots_purity_{ploidy}.svg",
    conda:
        "envs/R.yaml"
    script:
        "scripts/Compare_purity_measures.R"
