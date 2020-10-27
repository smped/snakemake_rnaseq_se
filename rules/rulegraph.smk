rule make_rulegraph:
    output:
        dot = "rules/rulegraph.dot",
        pdf = "rules/rulegraph.pdf"
    conda:
        "../envs/graphviz.yml"
    shell:
        """
        snakemake --rulegraph > {output.dot}
        dot -Tpdf {output.dot} > {output.pdf}
        """
