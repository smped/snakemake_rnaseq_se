rule make_rulegraph:
    output:
        dot = "output/rulegraph.dot",
        pdf = "output/rulegraph.pdf"
    conda:
        "../envs/graphviz.yml"
    shell:
        """
        snakemake --rulegraph > {output.dot}
        dot -Tpdf {output.dot} > {output.pdf}
        """
