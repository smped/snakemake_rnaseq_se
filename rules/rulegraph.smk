rule make_rulegraph:
    output:
        dot = "output/rulegraph.dot"
    shell:
        """
        snakemake --rulegraph > {output.dot}
        """

rule convert_rulegraph:
    input:
        dot = rules.make_rulegraph.output.dot
    output:
        pdf = "output/rulegraph.pdf"
    conda:
        "../envs/graphviz.yml"
    shell:
        """
        dot -Tpdf {input.dot} > {output.pdf}
        """
