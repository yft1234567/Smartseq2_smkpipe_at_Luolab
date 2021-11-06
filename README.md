# Smartseq2_smkpipe_at_Luolab
Snakemake pipeline for Upstream processing of Smartseq2 sequenced libraries

## Structure of the Repository

We follow a modified git repository structure snakemake pipelines as recommended by <snakmake docs>[https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html].

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
