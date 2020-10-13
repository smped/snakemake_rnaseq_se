rule get_genome:
    output:
        ref_path
    params:
        species = species,
        datatype = "dna",
        build = build,
        release = release
    log:
        "logs/refs/get_genome.log"
    wrapper:
        "0.66.0/bio/reference/ensembl-sequence"

rule get_annotation:
    output:
        gtf_path
    params:
        species = species,
        release = release,
        build = build,
        fmt = "gtf",
        flavor = "" # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        "logs/refs/get_annotation.log"
    wrapper:
        "0.66.0/bio/reference/ensembl-annotation"

rule star_index:
    input:
        fasta = rules.get_genome.output,
        gtf = rules.get_annotation.output
    output:
        directory(star_dir)
    conda:
        "../envs/star.yml"
    message:
        "Building STAR index"
    threads:
        16
    params:
        extra = config['star']['indexing_extra'],
        sjdbOverhang = config['star']['sjdbOverhang']
    log:
        "logs/refs/star_index.log"
    shell:
        """
        STAR \
          --runMode genomeGenerate \
          {params.extra} \
          --runThreadN {threads} \
          --genomeDir {output} \
          --genomeFastaFiles {input.fasta} \
          --sjdbOverhang {params.sjdbOverhang} \
          --sjdbGTFfile {input.gtf} &> {log}
        """
