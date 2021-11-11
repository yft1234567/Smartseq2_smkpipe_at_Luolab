# Step 1: Identify cell barcode whitelist (identify correct BC)
rule umi_tools_whitelist:
    input:
        read1="workflow/data/{user}/{project}/fastqs/{library}/{sample}_R2.fq.gz"
    output:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_whitelist.txt"
    log:
        "workflow/data/{user}/{project}/logs/{library}/{sample}_whitelist.log"
    conda:
        "../envs/umi_tools.yaml"
    threads:1
    shell:
        """
        umi_tools whitelist --bc-pattern=CCCCCCCCNNNNNNNN \
                            --log={log} \
                            --stdin {input.read1} \
                            --set-cell-number=100 \
                            --plot-prefix=workflow/data/{wildcards.user}/{wildcards.project}/alignment/{wildcards.library}/{wildcards.sample} \
                            --log2stderr > {output}
        """

# Step 2: Wash whitelist
rule wash_whitelist:
    input:
        whitelist="workflow/data/{user}/{project}/alignment/{library}/{sample}_whitelist.txt",
        ground_truth=config['ground_truth']
    output:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_whitelist_washed.txt"
    threads:1
    script:
        "../scripts/wash_whitelist.py"

# Step 3: Extract barcodes and UMIs and add to read names
rule umi_tools_extract:
    input:
        read1="workflow/data/{user}/{project}/fastqs/{library}/{sample}_R2.fq.gz",
        read2="workflow/data/{user}/{project}/fastqs/{library}/{sample}_R1.fq.gz",
        whitelist_washed="workflow/data/{user}/{project}/alignment/{library}/{sample}_whitelist_washed.txt"
    output:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_extracted.fq.gz"
    log:
        "workflow/data/{user}/{project}/logs/{library}/{sample}_extract.log"
    conda:
        "../envs/umi_tools.yaml"
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

# Step 4-0: Generate genome index
rule STAR_gen:
    input:
        genome_fa=config["genome_fa"],
        gtf=config["gtf_annotation"]
    output:
        directory(config["genome_index"])
    params:
        sjdbOverhang=config["sjdbOverhang"]
    conda:
        "../envs/star.yaml"
    threads:32
    shell:
        """
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome_fa} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang}
        """

# Step 4-1: Map reads
rule STAR:
    input:
        extracted_fq="workflow/data/{user}/{project}/alignment/{library}/{sample}_extracted.fq.gz",
        genomeDir=config["genome_index"]
    output:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_Aligned.sortedByCoord.out.bam"
    conda:
        "../envs/star.yaml"
    threads:16
    # STAR sometimes fails because of too many opened files (due to high thread count)
    # Lower thread here if necessary
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeLoad LoadAndKeep \
             --genomeDir {input.genomeDir} \
             --readFilesIn {input.extracted_fq} \
             --readFilesCommand zcat \
             --limitBAMsortRAM 20000000000 \
             --outFilterMultimapNmax 1 \
             --outFilterType BySJout \
             --outSAMstrandField intronMotif \
             --outFilterIntronMotifs RemoveNoncanonical \
             --outFilterMismatchNmax 6 \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outFileNamePrefix workflow/data/{wildcards.user}/{wildcards.project}/alignment/{wildcards.library}/{wildcards.sample}_
        """

# Step 4-2: Unload STAR genome index
rule STAR_unload:
    input:
        bams=get_files("STAR"),
        genomeDir=config["genome_index"]
    output:
        temp(touch("workflow/data/{user}/{project}/logs/STARunload.done")),
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --genomeLoad Remove \
             --genomeDir {input.genomeDir} \
             --outSAMmode None
        """

# Step 5-1: Assign reads to genes (featureCount)
rule featurecount:
    input:
        gtf=config["gtf_annotation"],
        bam="workflow/data/{user}/{project}/alignment/{library}/{sample}_Aligned.sortedByCoord.out.bam",
        dummy="workflow/data/{user}/{project}/logs/STARunload.done"
    output:
        assigned="workflow/data/{user}/{project}/alignment/{library}/{sample}_gene_assigned",
        bam_counted=temp("workflow/data/{user}/{project}/alignment/{library}/{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam")
    conda:
        "../envs/featurecount.yaml"
    threads:32
    shell:
        """
        featureCounts -s 1 \
                      -a {input.gtf} \
                      -o {output.assigned} \
                      -R BAM {input.bam} \
                      -T {threads}
        """

# Step 5-2: Assign reads to genes (sort bam files)
rule sambamba_sort:
    input:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam"
    output:
        temp("workflow/data/{user}/{project}/alignment/{library}/{sample}_assigned_sorted.bam")
    params:
        extra="-t "+str(config["thread_use"])+" -m 64G"
    threads: 32
    wrapper:
        "0.79.0/bio/sambamba/sort"

# Step 6: Count UMIs per gene per cell
rule umi_tools_count:
    input:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_assigned_sorted.bam"
    output:
        temp("workflow/data/{user}/{project}/alignment/{library}/{sample}_counts_raw.tsv.gz")
    conda:
        "../envs/umi_tools.yaml"
    threads:1
    shell:
        """
        umi_tools count --per-gene \
                        --gene-tag=XT \
                        --per-cell \
                        --stdin={input} \
                        --stdout={output}
        """

# Step 7: Append suffix to cells
rule append_sfx:
    input:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_counts_raw.tsv.gz"
    output:
        "workflow/data/{user}/{project}/alignment/{library}/{sample}_counts.tsv.gz"
    threads:1
    script:
        "../scripts/append_suffix.py"

# Step 8: Aggregate counts
rule aggr_counts:
    input:
        get_files("append_sfx")
    output:
        "workflow/data/{user}/{project}/outs/"+config["project"]+"_counts_all.tsv.gz"
    threads:1
    shell:
        "zcat {input} | gzip > {output}"

# # Step 9-1: Parse Seurat Object
# rule parse_seurat:
#     input:
#         "workflow/data/{user}/{project}/outs/"+config["project"]+"_counts_all.tsv.gz"
#     output:
#         "workflow/data/{user}/{project}/outs/"+config["project"]+"_seurat.rds"
#     threads:1
#     script:
#         "../scripts/parse_seurat.R"

# # Step 9-2: Parse AnnData Object
# rule parse_anndata:
#     input:
#         "workflow/data/{user}/{project}/outs/"+config["project"]+"_counts_all.tsv.gz"
#     output:
#         "workflow/data/{user}/{project}/outs/"+config["project"]+"_adata.h5ad"
#     threads:1
#     script:
#         "../scripts/parse_anndata.py"
