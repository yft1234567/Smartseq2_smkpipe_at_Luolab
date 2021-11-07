rule umi_tools_whitelist:
    input:
        read1="workflow/data/{user}/fastqs/{library}/{sample}_R1_001.fastq.gz"
    output:
        "workflow/data/{user}/alignment/{library}/{sample}_whitelist.txt"
    log:
        "workflow/data/{user}/logs/{library}/{sample}_whitelist.log"
    threads:1
    shell:
        """
        umi_tools whitelist --bc-pattern='CCCCCCCCNNNNNNNN' \
                            -L {log} \
                            --stdin {input.read1} \
                            --set-cell-number=100 \
                            --plot-prefix=workflow/data/zhongshilin/{wildcards.library}/outs/{wildcards.sample} \
                            --log2stderr > {output}
        """

rule wash_whitelist:
    input:
        whitelist="workflow/data/{user}/alignment/{library}/{sample}_whitelist.txt"
        ground_truth=""
    output:
        "workflow/data/{user}/alignment/{library}/{sample}_whitelist_washed.txt"
    threads:1
    script:
        "scripts/wash_whitelist.py"

rule umi_tools_extract:
    input:
        read1="workflow/data/{user}/fastqs/{library}/{sample}_R1_001.fastq.gz",
        whitelist="workflow/data/{user}/alignment/{library}/{sample}_whitelist_washed.txt"
    output:
        "workflow/data/{user}/{library}/fastqs/{sample}_R1_001_extracted.fastq.gz"
    log:
        "workflow/data/{user}/{library}/logs/{sample}_extract.log"
    threads:1
    shell:
        """
        umi_tools extract --extract-method=regex \
                          --bc-pattern='(?P<cell_1>.{{6}})CTTGTGGAAAGGACGAAACA{{s<=2}}(?P<cell_2>.{{6}})(?P<umi_1>.{{15}}).*' \
                          -L {log} \
                          --stdin {input.read1} \
                          --error-correct-cell \
                          --filter-cell-barcode \
                          --stdout {output} \
                          --whitelist={input.whitelist}
        """