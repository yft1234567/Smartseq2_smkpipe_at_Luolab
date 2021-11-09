# Smartseq2_smkpipe_at_Luolab
Snakemake pipeline for Upstream processing of Smartseq2 sequenced libraries

## Structure of the Repository

We follow a modified git repository structure snakemake pipelines as recommended by [snakmake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).

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
|   │   │   │   └── logs
|   │   │   └── Library2...
|   |   └── User2...
│   ├── rules
|   │   ├── module1.smk
|   │   └── module2.smk
│   ├── envs
|   │   ├── tool1.yaml
|   │   └── tool2.yaml
│   ├── scripts
|   │   ├── script1.py
|   │   └── script2.R
│   ├── notebooks
|   │   ├── notebook1.py.ipynb
|   │   └── notebook2.r.ipynb
│   ├── report
|   │   ├── plot1.rst
|   │   └── plot2.rst
|   └── Snakefile
├── config
│   ├── config.yaml
│   └── some-sheet.tsv
```

## Pipeline details

We thankfully adopt a comprehensive pipeline described in the [umi_tools documentation](https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html). We made some custom modification to the pipeline for the purpose of speed boost by exploiting multi-threading features and parallel processing.

```
#! /bin/env bash
# Step 1: Identify correct cell barcodes (umi_tools whitelist)
umi_tools whitelist --bc-pattern=CCCCCCCCNNNNNNNN \
		    --log=whitelist.log \
		    --stdin=SAMPLE_R2.fq.gz
		    --set-cell-number=100 \
		    --plot-prefix=SAMPLE
		    --log2stderr > whitelist.txt

# Step 2: Wash whitelist
(run scripts/wash_whitelist.py) 

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

# Step 4: Map reads
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

# Step 5-1: Assign reads to genes
featureCounts -s 1 \
	      -a annotation.gtf \
	      -o SAMPLE_gene_assigned \
	      -R BAM SAMPLE_Aligned.sortedByCoord.out.bam \
	      -T 32

# Step 5-2: Sort and index BAM file
sambamba sort -t 32 -m 64G 

# Step 6: Count UMIs per gene per cell
umi_tools count --per-gene \
		--gene-tag=XT \
		--per-cell \
		--stdin=SAMPLE_assigned_sorted.bam \
		--stdout=SAMPLE_counts.tsv.gz

# Step 7: Unload STAR genome
STAR --genomeLoad Remove \
     --genomeDir /path/to/genome/index

```

