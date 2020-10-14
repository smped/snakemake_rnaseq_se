rule make_rproj:
    output: os.getcwd() + ".Rproj"
    threads: 1
    shell:
        """
        if [[ ! -f {output} ]]; then
          echo -e "Version: 1.0\n" > {output}
          echo -e "RestoreWorkspace: Default\nSaveWorkspace: Default\nAlwaysSaveHistory: Default\n" >> {output}
          echo -e "EnableCodeIndexing: Yes\nUseSpacesForTab: Yes\nNumSpacesForTab: 2\nEncoding: UTF-8\n" >> {output}
          echo -e "RnwWeave: knitr\nLaTeX: pdfLaTeX\n" >> {output}
          echo -e "AutoAppendNewline: Yes\nStripTrailingWhitespace: Yes" >> {output}
        fi
        """

rule build_wflow_description:
    input:
        dot = rules.make_rulegraph.output.dot,
        rmd = "analysis/description.Rmd",
        rproj = rules.make_rproj.output
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
        rmd = "analysis/qc_raw.Rmd",
        rproj = rules.make_rproj.output
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
        rmd = "analysis/qc_trimmed.Rmd",
        rproj = rules.make_rproj.output
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
        rmd = "analysis/qc_aligned.Rmd",
        rproj = rules.make_rproj.output
    output:
        html = "docs/qc_aligned.html",
        rds = "output/genesGR.rds"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/qc_aligned.log"
    threads: 1
    shell:
       """
       R -e "workflowr::wflow_build('{input.rmd}')" 2>&1 > {log}
       """

rule build_dge_analysis:
    input:
        rds = rules.build_qc_aligned.output.rds,
        rmd = "analysis/dge_analysis.Rmd",
        rproj = rules.make_rproj.output
    output:
        html = "docs/dge_analysis.html"
    conda:
        "../envs/workflowr.yml"
    log:
        "logs/workflowr/dge_analysis.log"
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
        aligned = rules.build_qc_aligned.output.html,
        dge = rules.build_dge_analysis.output.html,
        rproj = rules.make_rproj.output
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
