configfile: "config.yaml"
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]

(subtypes,projects,hash_samples,files,) = glob_wildcards(data_dir + "{subtype}/{project}/harmonized/Copy_Number_Variation/Masked_Copy_Number_Segment/{hash_sample}/{file}")

#-------------------------------------------------------------------------------------------------------------------
# 0.3 specify target rules
rule all:
    input:
        

#++++++++++++++++++++++++++++++++++++++++++++++++ 1 DOWNLOAD DATA  +++++++++++++++++++++++++++++++++++++++++++++++++
