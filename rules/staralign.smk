rule star_se:
    input:
        fq1 = "data/trimmed/fastq/{sample}" + ext,
        index = rules.star_index.output
    output:
        "data/aligned/bam/{sample}/Aligned.sortedByCoord.out.bam"
    conda:
        "../envs/star.yml"
    log:
        "logs/star/{sample}.log"
    params:
        extra = "--outSAMtype BAM SortedByCoordinate"
    threads: 8
    script:
        "../scripts/star_alignment.py"

rule index_bam:
    input:
        rules.star_se.output
    output:
        "data/aligned/bam/{sample}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/samtools.yml"
    threads: 1
    shell:
        """
        samtools index {input} {output}
        """
