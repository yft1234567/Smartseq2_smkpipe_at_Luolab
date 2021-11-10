import pandas as pd

def append_suffix(file_path,
                  sample):
    counts = pd.read_csv(file_path, sep="\t")
    counts['cell'] = counts['cell'] + '_' + sample
    return counts

counts = append_suffix(file_path = snakemake.input[0],
                       sample = snakemake.wildcards.sample)
counts.to_csv(snakemake.output[0], index=False, sep="\t", compression="gzip")