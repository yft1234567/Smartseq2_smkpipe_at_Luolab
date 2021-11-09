import pandas as pd
# import glob
# import os

# from snakemake.utils import validate

# validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], dtype={'User': str, 'Project': str, 'Sample': str})
    .set_index("Sample", drop=False)
    .sort_index()
)

def parse_suffix(rule):
    # Given a rule name, return the suffix of its corresponding output file
    match rule:
        case 'umi_tools_whitelist':
            return 'whitelist.txt'
        case 'wash_whitelist':
            return 'whitelist_washed.txt'
        case 'umi_tools_extract':
            return 'extracted.fq.gz'
        case 'STAR':
            return 'Aligned.sortedByCoord.out.bam'
        case 'featurecount':
            return 'Aligned.sortedByCoord.out.bam.featureCounts.bam'
        case 'sambamba_sort':
            return 'assigned_sorted.bam'
        case 'umi_tools_count':
            return 'counts.tsv.gz'


def get_files(rule):

    files = expand(
        '_'.join(["workflow/data/{user}/{project}/alignment/{sample}/{sample}_", parse_suffix(rule)])
        zip, user=samples.User.to_list(), project=samples.Project.to_list(), sample=samples.Sample.to_list()
        )
    return files