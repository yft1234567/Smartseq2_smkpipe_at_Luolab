# Step 1: Identify cell barcode whitelist (identify correct BC)
rule umi_tools_whitelist:
    input:
        read1="workflow/data/{user}/fastqs/{library}/{sample}_sub_R2.fq.gz"
    output:
        "workflow/data/{user}/alignment/{library}/{sample}_whitelist.txt"
    log:
        "workflow/data/{user}/logs/{library}/{sample}_whitelist.log"
    threads:1
    shell:
        """
        umi_tools whitelist --bc-pattern=CCCCCCCCNNNNNNNN \
                            -L {log} \
                            --stdin {input.read1} \
                            --set-cell-number=100 \
                            --plot-prefix=workflow/data/{wildcards.user}/alignment/{wildcards.library}/{wildcards.sample} \
                            --log2stderr > {output}
        """

# Step 2: Wash whitelist
rule wash_whitelist:
    input:
        whitelist="workflow/data/{user}/alignment/{library}/{sample}_whitelist.txt",
        ground_truth=config['ground_truth']
    output:
        "workflow/data/{user}/alignment/{library}/{sample}_whitelist_washed.txt"
    threads:1
    script:
        "../scripts/wash_whitelist.py"

# Step 3: Extract batcodes and UMIs and add to read names
rule umi_tools_extract:
    input:
        read1="workflow/data/{user}/fastqs/{library}/{sample}_sub_R2.fq.gz"
        read2="workflow/data/{user}/fastqs/{library}/{sample}_sub_R1.fq.gz"
        whitelist_washed="workflow/data/{user}/alignment/{library}/{sample}_whitelist_washed.txt"
    output:
        "workflow/data/{user}/alignment/{library}/{sample}_sub_extracted.txt"
    log:
        "workflow/data/{user}/logs/{library}/{sample}_extract.log"
    threads:1
    shell:
        """
        umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                          --log={log} \
                          --stdin {input.read1} \
                          --read2-in {input.read2} \
                          --stdout {output} \
                          --read2-stdout \
                          --filter-cell-barcode \
                          --error-correct-cell \
                          --whitelist={input.whitelist_washed}
        """

# # Step 4: Map reads
# rule STAR:
#     input:
#         extracted_fq="workflow/data/{user}/alignment/{library}/{sample}_sub_extracted.txt"
#         genomeDir="/mnt/ZA1BT1ER/raywang/ensembl/mus_musculus/STAR_INDEX_GRCm38_102_mScarlet/"
#     output:
#         "workflow/data/{user}/alignment/{library}/{sample}_"
#     log:

#     threads:32
#     shell:
#         """
#         STAR --runThreadN 32 \
#              --genomeDir {input.genomeDir} \
#              --readFilesIn {input.extracted_fq} \
#              --readFilesCommand zcat \
#              --outFilterMultimapNmax 1 \
#              --outFilterType BySJout \
#              --outSAMstrandField intronMotif \
#              --outFilterIntronMotifs RemoveNoncanonical \
#              --outFilterMismatchNmax 6 \
#              --outSAMtype BAM SortedByCoordinate \
#              --outFileNamePrefix {output} \
#              --outSAMunmapped Within
#         """

# # Step 5: Assign reads to genes (featureCount)
# rule featurecount:
#     input:
    