rule umi_tools_whitelist:
    input:
        "workflow/data/{user}/fastqs/{library}/{sample}_R1_001.fastq.gz"
    output:
        "workflow/data/{user}/alignment/{library}/{sample}_whitelist.txt"
    log:
        "workflow/data/{user}/logs/{library}/{sample}_whitelist.log"
    threads:1
    shell:
        """
        umi_tools whitelist --extract-method=regex \
                            --bc-pattern='(?P<cell_1>.{{6}})CTTGTGGAAAGGACGAAACA{{s<=2}}(?P<cell_2>.{{6}})(?P<umi_1>.{{15}}).*' \
                            -L {log} \
                            --stdin {input} \
                            --set-cell-number=800 \
                            --plot-prefix=workflow/data/zhongshilin/{wildcards.library}/outs/{wildcards.sample} \
                            --log2stderr > {output}
        """

rule umi_tools_extract:
    input:
        fq="workflow/data/zhongshilin/{library}/fastqs/{sample}_R1_001.fastq.gz",
        whitelist="workflow/data/zhongshilin/{library}/outs/{sample}_whitelist_washed.txt"
    output:
        "workflow/data/zhongshilin/{library}/fastqs/{sample}_R1_001_extracted.fastq.gz"
    log:
        "workflow/data/zhongshilin/{library}/logs/{sample}_extract.log"
    threads:1
    shell:
        """
        umi_tools extract --extract-method=regex \
                          --bc-pattern='(?P<cell_1>.{{6}})CTTGTGGAAAGGACGAAACA{{s<=2}}(?P<cell_2>.{{6}})(?P<umi_1>.{{15}}).*' \
                          -L {log} \
                          --stdin {input.fq} \
                          --error-correct-cell \
                          --filter-cell-barcode \
                          --stdout {output} \
                          --whitelist={input.whitelist}
        """