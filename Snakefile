import pandas as pd
import os.path

configfile: "config/config.yml"

# Config information for the reference
species = config['ref']['species']
build = config['ref']['build']
release = str(config['ref']['release'])
seqtype = config['ref']['seqtype']

# Define the path for the reference
ens_full = "ensembl-release-" + release
ref_root = os.path.join(config['ref']['root'], ens_full, species)
ref_fa = species.capitalize() + "." + build + "." + "dna." + seqtype + ".fa"
ref_path = os.path.join(ref_root, 'dna', ref_fa)

# Define the path for the GTF annotation
# By default, the wrapper will extract the file
gtf = species.capitalize() + "." + build + "." + release + ".gtf"
gtf_path = os.path.join(ref_root, gtf)
sj_oh = str(config['star']['sjdbOverhang'])
star_dir = os.path.join(ref_root, 'dna', 'star_sjOverhang' + sj_oh)

# Samples
ext = config['ext']
tag = config['tag']
samples = pd.read_table(config["samples"])
counts_file = "data/aligned/counts/counts.out"

# Define outputs
ALL_RULEGRAPH = expand(["rules/rulegraph.{suffix}"],
                       suffix = ['dot', 'pdf'])
ALL_REFS = [ref_path, gtf_path, star_dir]
ALL_FQC = expand(["data/{step}/FastQC/{sample}{tag}_fastqc.{suffix}"],
                 suffix = ['zip', 'html'],
                 tag = [tag],
                 sample = samples['sample'],
                 step = ['raw', 'trimmed'])
ALL_TRIMMED = expand(["data/trimmed/fastq/{sample}" + tag + ext],
                      sample = samples['sample'])
ALL_ALN = expand(["data/aligned/bam/{sample}{tag}/Aligned.sortedByCoord.out.{suffix}"],
                 suffix = ['bam.bai'],
                 tag = [tag],
                 sample = samples['sample'])
ALL_WORKFLOWR = expand(["docs/{step}.html"],
                      step = ['description', 'index', 'qc_raw'])
ALL_OUTPUTS = []
ALL_OUTPUTS.extend(ALL_RULEGRAPH)
ALL_OUTPUTS.extend(ALL_REFS)
ALL_OUTPUTS.extend(ALL_TRIMMED)
ALL_OUTPUTS.extend(ALL_FQC)
ALL_OUTPUTS.extend(ALL_ALN)
ALL_OUTPUTS.extend(ALL_WORKFLOWR)
ALL_OUTPUTS.extend([counts_file])

# And the rules
rule all:
    input:
        ALL_OUTPUTS

include: "rules/refs.smk"
include: "rules/rulegraph.smk"
include: "rules/trimming.smk"
include: "rules/qc.smk"
include: "rules/staralign.smk"
include: "rules/count.smk"
include: "rules/workflowr.smk"
