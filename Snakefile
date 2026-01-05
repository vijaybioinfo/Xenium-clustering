"""
### Snakefile for Spatial transcriptomics
"""
import pandas as pd
import os
import re

conf_file_loc = "config.yaml"

configfile:conf_file_loc

# import yaml
# with open('/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/spatial_pilot_5k/scripts/devel/config.yaml', 'r') as file:
#     config = yaml.safe_load(file)

# file_path = config["input_expression"] +  "transcripts.csv.gz"

# if os.path.exists(file_path):
#     print(f"Transcripts file '{file_path}' exists.")
# else:
#     print(f"Transcripts file don't exist. This will take som extra time '{file_path}'")

get_transcripts_csv_script=config["pipeline"] + "/R/get_transcripts_csv.R"
init_script=config["pipeline"] + "/R/init.R"
markers_script=config["pipeline"] + "/R/markers.R"
compo_script=config["pipeline"] + "/R/report_components.R"

MEANS = config['variable_features']['mean.cutoff'][0]
PERCENTAGES = config['variable_features']['percent']
COMPONENTES = config['dim_reduction']['base']['chosen_comp']
RESOLUTIONS = config['resolution']
EXECS = config['exec']

rule all:
    input:
        ".transcripts_done.txt",
         expand(".object_init_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
         expand(".object_metadata_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
         expand(".object_reductions_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
         expand(".object_graphs_mean{mean}_pct{percentage}_pc{component}.rds", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES),
        expand(".markers_mean{mean}_pct{percentage}_pc{component}_res{resolution}.txt", mean = MEANS ,percentage = PERCENTAGES, component = COMPONENTES, resolution = RESOLUTIONS),
        expand(".report_mean{mean}_pct{percentage}_pc{component}.txt", mean = MEANS, percentage = PERCENTAGES, component = COMPONENTES)

 ### --------------- Cellranger parquet file prepocessing --------------- ###
rule get_transcripts_csv:
    input:
        conf_file_loc
    output:
        ".transcripts_done.txt"
    shell: 
        '{EXECS} {get_transcripts_csv_script} --yaml {conf_file_loc} -v TRUE'

 ### --------------- Seurat Normalization and HVG selection --------------- ###
rule init_object:
    input:
        ".transcripts_done.txt"
    output:
        ".object_init_mean{mean}_pct{percentage}_pc{component}.rds",
        ".object_metadata_mean{mean}_pct{percentage}_pc{component}.rds",
        ".object_reductions_mean{mean}_pct{percentage}_pc{component}.rds",
        ".object_graphs_mean{mean}_pct{percentage}_pc{component}.rds"
    params:
        component = config['dim_reduction']['base']['n_comp']
    shell:
        "{EXECS} {init_script} -y {conf_file_loc} --percent {wildcards.percentage} --n_comp {params.component}  --chosen_comp {wildcards.component} --prefix init_mean{wildcards.mean}_pct{wildcards.percentage}_pc{wildcards.component}"

rule markers:
    input:
        ".object_init_mean{mean}_pct{percentage}_pc{component}.rds"
    output:
        ".markers_mean{mean}_pct{percentage}_pc{component}_res{resolution}.txt"
    message: " --- Branch resolution for marker calculation --- "
    shell:
        "{EXECS} {markers_script} -y {conf_file_loc} --init {input} --percent {wildcards.percentage} --chosen_comp {wildcards.component} --resolution {wildcards.resolution} --prefix init_mean{wildcards.mean}_pct{wildcards.percentage}_pc{wildcards.component}_res{wildcards.resolution}"

rule report_components:
    input:
        ".object_init_mean{mean}_pct{percentage}_pc{component}.rds"
    output:
        ".report_mean{mean}_pct{percentage}_pc{component}.txt"
    message: " --- Creating report: components  ---"
    shell:
        "{EXECS} {compo_script} -y {conf_file_loc} --init {input} --percent {wildcards.percentage} --chosen_comp {wildcards.component} --prefix init_mean{wildcards.mean}_pct{wildcards.percentage}_pc{wildcards.component}"