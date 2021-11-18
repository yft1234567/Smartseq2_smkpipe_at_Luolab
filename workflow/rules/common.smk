import os
import pandas as pd
# import glob


# from snakemake.utils import validate

# validate(config, schema="../schemas/config.schema.yaml")
## Debug
# import yaml
# with open ('config/config.yaml') as f:
#     config = yaml.safe_load(f)

samples = (
    pd.read_csv(config["samples"], dtype={'User': str, 'Project': str, 'Library': str, 'Sample': str})
    .set_index("Sample", drop=False)
    .sort_index()
)

samples.Sample.fillna(samples.Library, inplace=True)

def parse_suffix(rule):
    # Given a rule name, return the suffix of its corresponding output file
    if rule == 'umi_tools_whitelist':
        return 'whitelist.txt'
    elif rule == 'wash_whitelist':
        return 'whitelist_washed.txt'
    elif rule == 'umi_tools_extract':
        return 'extracted.fq.gz'
    elif rule == 'STAR':
        return 'Aligned.sortedByCoord.out.bam'
    elif rule == 'featurecount':
        return 'Aligned.sortedByCoord.out.bam.featureCounts.bam'
    elif rule == 'sambamba_sort':
        return 'assigned_sorted.bam'
    elif rule == 'umi_tools_count':
        return 'counts_raw.tsv.gz'
    elif rule == "append_sfx":
        return 'counts.tsv.gz'

def get_files(rule):
    files = expand(
        '_'.join(["workflow/data/{user}/{project}/alignments/{library}/{sample}", parse_suffix(rule)]),
        zip, user=samples.User.to_list(), project=samples.Project.to_list(), library=samples.Library.to_list(), sample=samples.Sample.to_list()
        )
    return files

def get_aggr_output():
    # return "workflow/data/LinRui/AAV_Microglia/outs/exp_mat.tsv.gz"
    return list(set(expand("workflow/data/{user}/{project}/outs/{project}_counts_all.tsv.gz", 
    zip, user=samples.User.to_list(), project=samples.Project.to_list())))