# Smartseq2 Snakemake Pipeline at Luolab
Snakemake pipeline for Upstream processing of Smartseq2 sequenced libraries - **from fastq.gz to count table**.

## Requirements

If you choose to use `Integrated Package Management` feature of Snakemake, you can skip this section. The pipeline will automatically builds most of the required tools listed in the `workflow/envs/...yaml` files when you activate the `--use-conda` flag. But note you should still install `python` and `pandas` as a minimum requirement.

Otherwise, you may install the following packages into your environment before running the pipeline. Go to the manual pages of these packages and follow the installation instructions to build them properly.

```
python (>=3.1.0)
pandas (>=1.1.0)
umi_tools (>=1.1.2)
STAR (>=2.5.4)
subread (>=1.6.3)
```

## Usage

Once your environment is built, edit `sample_table.csv` and `config.yaml` files in the `config/` directory. Depending on your system specifications, you may also need to manually edit some parameters in the workflow (i.e. `pipeline.smk`).

Activate a shell prompt, change to main directory of the pipeline (i.e. `Smartseq2_smkpipe_at_Luolab/`) and type

```
snakemake --cores XX
```

or

```
snakemake --cores XX --use-conda
```

if you wish an automatic package management.

Make sure to use the correct genome index, barcode-ground-truth list and gtf annotations, as specified in the `config/config.yaml`. Also make sure to construct the correct file structure under the `data/` directory. 

For proper file system management, we recommend to put raw fastq files (`fastqs/`), intermediate files (`alignment/`), and aggregated final results (`aggr_outs/`) in dedicated locations in your hardware, and construct symbolic links between your `data/` directory and the file storage systems, as shown in the project structure below.

## Structure of the Repository

We follow a git repository structure snakemake pipelines as recommended by [snakmake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) with minor modifications.


```
├── .gitignore
├── README.md
├── LICENSE.md
├── workflow
|   ├── data
|   |   ├── User1
|   │   │   ├── Library1
|   │   │   │   ├── fastqs
|   │   │   │   ├── alignment
|   │   │   │   ├── aggr_outs
|   │   │   │   └── logs
|   │   │   └── Library2
|   |   └── User2
│   ├── rules
|   │   ├── common.smk
|   │   └── pipeline.smk
│   ├── envs
|   │   ├── tool1.yaml
|   │   └── tool2.yaml
│   ├── scripts
|   │   ├── wash_whitelist.py
|   │   └── script2.R
│   ├── notebooks
|   │   ├── notebook1.py.ipynb
|   │   └── notebook2.r.ipynb
│   ├── report
|   │   ├── plot1.rst
|   │   └── plot2.rst
|   └── Snakefile
└── config
    ├── config.yaml
    └── sample_table.tsv
```

## Pipeline details

We adopt a comprehensive pipeline described in the [umi_tools documentation](https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html). Some custom modification were made for speed boosting (i.e. multi-threading, parallel processing, STAR "shared memory" modules).

```
#! /bin/env bash
# Step 1: Identify correct cell barcodes (umi_tools whitelist)
umi_tools whitelist --bc-pattern=CCCCCCCCNNNNNNNN \
		    --log=whitelist.log \
		    --stdin=SAMPLE_R2.fq.gz \  ## R1 and R2 are swapped b.c. we customized our protocol
		    --set-cell-number=100 \
		    --plot-prefix=SAMPLE \
		    --log2stderr > whitelist.txt

# Step 2: Wash whitelist
(impletemded in scripts/wash_whitelist.py) 

# Step 3: Extract barcodes and UMIs and add to read names
umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
		  --log=extract.log \
		  --stdin=SAMPLE_R2.fq.gz \
		  --read2-in=SAMPLE_R1.fq.gz \
		  --stdout=SAMPLE_extracted.fq.gz \
		  --read2-stdout \
		  --filter-cell-barcode \
		  --error-correct-cell \
		  --whitelist=whitelist_washed.txt

# Step 4-1: Map reads
STAR --runThreadN 32 \
     --genomeLoad LoadAndKeep \
     --genomeDir /path/to/genome/index \
     --readFilesIn SAMPLE_extracted.fq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outFilterType BySJout \
     --outSAMstrandField intronMotif \
     --outFilterIntronMotifs RemoveNoncanonical \
     --outFilterMismatchNmax 6 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFileNamePrefix SAMPLE_

# Step 4-2: Unload STAR genome
STAR --genomeLoad Remove \
     --genomeDir /path/to/genome/index

# Step 5-1: Assign reads to genes
featureCounts -s 1 \
	      -a annotation.gtf \
	      -o SAMPLE_gene_assigned \
	      -R BAM SAMPLE_Aligned.sortedByCoord.out.bam \
	      -T 32

# Step 5-2: Sort and index BAM file
sambamba sort -t 32 -m 64G \
              -o SAMPLE_assigned_sorted.bam \
			  SAMPLE_Aligned.sortedByCoord.out.bam.featureCounts.bam 

# Step 6: Count UMIs per gene per cell
umi_tools count --per-gene \
		--gene-tag=XT \
		--per-cell \
		--stdin=SAMPLE_assigned_sorted.bam \
		--stdout=SAMPLE_counts.tsv.gz

# Step 7: Append suffix to cells
(implemented by scripts/append_suffix.py)

# Step 8: Aggregate counts
zcat SAMPLE1_counts.tsv.gz SAMPLE2_counts.tsv.gz ... | gzip > outs/counts_all.tsv.gz

```

## Notes

1. When executing the pipeline, depending on the system's hardware one may receive error: 
   > exiting because of OUTPUT FILE error: could not create output file SampleName-STARtmp//BAMsort/... SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.

   This is because STAR creates lots of temporary files during running and could reach file creation limit by the system. A simple workaround is to set the `threads:` in `pipeline.smk` `rule STAR: (# Step 4-1: Map reads)` to lower values, i.e. 32->16.

2. Within a project, we recommended that you fix your STAR version. This is because genomes generated by different versions of STARs are incompatible with each other. For example, an alignment job run by STAR v2.7.9 will not accept a genome index generated by STAR v2.6.1.