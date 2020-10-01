rule adapter_removal:
    input:
        r1 = "data/raw/fastq/{sample}" + ext
    output:
        r1 = "data/trimmed/fastq/{sample}" + ext,
        log = "data/trimmed/logs/{sample}.settings"
    conda:
        "../envs/adapterremoval.yml"
    params:
        adapter1 = config['trimming']['adapter1'],
        minlength = config['trimming']['minlength'],
        minqual = config['trimming']['minqual'],
        maxns = config['trimming']['maxns'],
        extra = config['trimming']['extra']
    threads: 1
    log:
        "logs/adapterremoval/{sample}.log"
    shell:
        """
        AdapterRemoval \
            --adapter1 {params.adapter1} \
            --file1 {input.r1} \
            --threads {threads} \
            {params.extra} \
            --maxns {params.maxns} \
            --minquality {params.minqual} \
            --minlength {params.minlength} \
            --output1 {output.r1} \
            --discarded /dev/null \
            --settings {output.log} &> {log}
        """
