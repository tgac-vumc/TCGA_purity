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
        expand("../data/copynumber/segments/{sample}_segments.txt", sample = Samples)
        

#++++++++++++++++++++++++++++++++++++++++++++++++ 1 PREPARE DATA  +++++++++++++++++++++++++++++++++++++++++++++++++

rule Migrate_files:
    input:
        paths
    output:
        dynamic("../data/copynumber/segments/{sample}_segments.txt")
    run:
        Migrate_files(Sample_dict, input, data_dir)
    
