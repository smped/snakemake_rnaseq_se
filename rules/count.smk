rule count:
    input:
        bams = expand(["data/aligned/bam/{sample}/Aligned.sortedByCoord.out.bam"],
               sample = samples['sample']),
        gtf = rules.get_annotation.output
    output:
        counts_file
    conda:
        "../envs/subread.yml"
    threads: 4
    params:
        fracOverlap = config['featureCounts']['fracOverlap'],
        q = config['featureCounts']['minQual'],
        s = config['featureCounts']['strandedness'],
        extra = config['featureCounts']['extra']
    shell:
       """
       featureCounts \
         {params.extra} \
         -Q {params.q} \
         -s {params.s} \
         --fracOverlap {params.fracOverlap} \
         -T {threads} \
         -a {input.gtf} \
         -o {output} \
         {input.bams}
       """
