rule make_rulegraph:
    output:
        dot = "rules/rulegraph.dot"
    shell:
        """
        snakemake --rulegraph > {output.dot}
        """

rule convert_rulegraph:
    input:
        dot = rules.make_rulegraph.output.dot
    output:
        pdf = "rules/rulegraph.pdf"
    conda:
        "../envs/graphviz.yml"
    shell:
        """
        dot -Tpdf {input.dot} > {output.pdf}
        """
