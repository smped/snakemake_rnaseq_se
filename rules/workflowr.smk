rule build_wflow_description:
    input:
        dot = rules.make_rulegraph.output.dot,
        rmd = "analysis/description.Rmd"
    output:
        html = "docs/description.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/description.log"
    threads: 1
    shell:
       """
       R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
       """

rule build_qc_raw:
    input:
        fqc = expand(["data/raw/FastQC/{sample}{tag}_fastqc.zip"],
                     tag = [tag], sample = samples['sample']),
        rmd = "analysis/qc_raw.Rmd"
    output:
        html = "docs/qc_raw.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/qc_raw.log"
    threads: 1
    shell:
       """
       R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
       """

rule build_qc_trimmed:
    input:
        fqc = expand(["data/trimmed/FastQC/{sample}{tag}_fastqc.zip"],
                     tag = [tag], sample = samples['sample']),
        rmd = "analysis/qc_trimmed.Rmd"
    output:
        html = "docs/qc_trimmed.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/qc_trimmed.log"
    threads: 1
    shell:
       """
       R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
       """

rule build_qc_aligned:
    input:
        counts = rules.count.output,
        aln_logs = expand(["data/aligned/bam/{sample}{tag}/Log.final.out"],
                          sample = samples['sample'], tag = [tag]),
        rmd = "analysis/qc_aligned.Rmd"
    output:
        html = "docs/qc_aligned.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/qc_aligned.log"
    threads: 1
    shell:
       """
       R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
       """

rule build_wflow_site_index:
    input:
        rmd = "analysis/index.Rmd",
        desc = rules.build_wflow_description.output.html,
        raw = rules.build_qc_raw.output.html,
        trimmed = rules.build_qc_trimmed.output.html,
        aligned = rules.build_qc_aligned.output.html
    output:
        html = "docs/index.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/index.log"
    threads: 1
    shell:
       """
       R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
       """
