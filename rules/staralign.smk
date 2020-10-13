rule star_se:
    input:
        fq1 = "data/trimmed/fastq/{sample}" + ext,
        index = rules.star_index.output
    output:
        bam = "data/aligned/bam/{sample}/Aligned.sortedByCoord.out.bam",
        log = "data/aligned/bam/{sample}/Log.final.out"
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
        rules.star_se.output.bam
    output:
        "data/aligned/bam/{sample}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/samtools.yml"
    threads: 1
    shell:
        """
        samtools index {input} {output}
        """
