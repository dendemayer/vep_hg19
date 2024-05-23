configfile: 'config.yaml'

rule all:
    input:
        #expand("annotated/{annotater}/{chr}.imputed.poly_subset.vcf.gz", chr=config['chromosomes'], annotater=['vep', 'snpeff']),
        "annotated/multiqc_report.html",

rule download_md5sums:
    output:
        md5sums="ref_seq/MD5SUMS",
    log:
        "logs/download_md5sums.log",
    benchmark:
        "benchmarks/download_md5sums.txt"
    conda: 
        "envs/wget.yaml",
    shell:
        "wget --no-clobber --directory-prefix=ref_seq https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/MD5SUMS"

rule download_refseq:
    '''
        the wget env yaml needs beside of wget also gzip to make sure thate the
        -k option is available
    '''
    input: 
        md5sums="ref_seq/MD5SUMS"
    output:
        gz="ref_seq/{ref_file}.gz",
        ref_file="ref_seq/{ref_file}",
    log:
        "logs/download_ref_seq/{ref_file}_all.log",
    benchmark:
        "benchmarks/download_ref_seq/{ref_file}.txt"
    wildcard_constraints:
        ref_file=".+fa|.+annotation.gtf|.+annotation.gff3"
    conda: 
        "envs/wget.yaml",
    shell:
        "scripts/download_refs.sh {wildcards.ref_file}.gz 2>> {log} && "
        "if [[ {output.gz} == *.gz ]]; then gzip -dkf {output.gz}; fi"
        #"wget --directory-prefix ref_seq {params.ref_link} 2> {log} && "
        #"gzip -d {output.gz} && bgzip {output.fa} 2>> {log}"

rule sort_gff_file:
    input:
        "ref_seq/{annot_file}.gff3"
    output:
        sort="ref_seq/{annot_file}_sorted.gff3",
        gz="ref_seq/{annot_file}_sorted.gff3.gz",
        tbi="ref_seq/{annot_file}_sorted.gff3.gz.tbi"
    conda:
        "envs/gff3sort.yaml"
    benchmark:
        "benchmarks/gff3_sort/{annot_file}.gff3.all"
    log:
        "logs/gff3_sort/{annot_file}.gff3.log"
    threads:
        8
    shell:
        "gff3sort.pl {input} > {output.sort} && "
        "bgzip -kf -@ {threads} {output.sort} && tabix -p gff {output.gz}"

rule samtools_index_ref:
    input:
        "ref_seq/{ref_file}.fa",
    output:
        "ref_seq/{ref_file}.fa.fai",
    log:
        "logs/samtools_index_ref_seq/{ref_file}_all.log",
    benchmark:
        "benchmarks/samtools_index_ref_seq/{ref_file}.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input} 2> {log}"

rule vep_annotate_variants_wrapper:
    input:
        calls="vcfs/{chr}.imputed.poly_subset.vcf.gz",  # .vcf, .vcf.gz or .bcf
        #cache="resources/vep/cache",  # can be omitted if fasta and gff are specified
        plugins="resources/VEP_plugins",
        # optionally add reference genome fasta
        fasta=expand("ref_seq/{ref_file}.fa", ref_file=config["ref_file"]),
         #fasta="genome.fasta",
        # fai="genome.fasta.fai", # fasta index
        fai=expand("ref_seq/{ref_file}.fa.fai", ref_file=config["ref_file"]),
        # gff from:
        #https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
        #gff=expand("ref_seq/{ref_file}.gff", ref_file=config["ref_file"]),
        gff=expand("ref_seq/{annotation_file}_sorted.gff3.gz", annotation_file=config["annotation_file"]),
        #gff="annotation.gff",
        # csi="annotation.gff.csi", # tabix index
        # add mandatory aux-files required by some plugins if not present in the VEP plugin directory specified above.
        # aux files must be defined as following: "<plugin> = /path/to/file" where plugin must be in lowercase
        # revel = path/to/revel_scores.tsv.gz
    output:
        calls="annotated/vep/{chr}.imputed.poly_subset.vcf.gz",  # .vcf, .vcf.gz or .bcf stats="vep/variants.html",
        #--stats_file [SAMPLE_NAME].vep.html (without the vep or summary suffix, MultiQC will ignore the HTML files)
        stats='annotated/vep/{chr}.imputed.poly_subset.vep.html',
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra="--everything",  # optional: extra arguments
    benchmark: 
        "benchmarks/vep/{chr}_benchmark.log"
    log:
        "logs/vep/{chr}.log",
    threads: 40
    wrapper:
        "v3.10.2/bio/vep/annotate"

rule snpeff_download:
    output:
        # wildcard {reference} may be anything listed in `snpeff databases`
        directory("resources/snpeff/{reference}")
    log:
        "logs/snpeff/download/{reference}.log"
    params:
        reference="{reference}"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "v3.10.2/bio/snpeff/download"

rule snpeff_wrapper:
    input:
        calls="vcfs/{chr}.imputed.poly_subset.vcf.gz",  # .vcf, .vcf.gz or .bcf
        db="resources/snpeff/GRCh37.75" # path to reference db downloaded with the snpeff download wrapper
    output:
        stats="annotated/snpeff/{chr}.imputed.poly_subset.html",  # summary statistics (in HTML), optional
        calls="annotated/snpeff/{chr}.imputed.poly_subset.vcf.gz",   # annotated calls (vcf, bcf, or vcf.gz)
        csvstats="annotated/snpeff/{chr}.imputed.poly_subset.csv", # summary statistics in CSV, optional
    log:
        "logs/snpeff/{chr}.log"
    resources:
        java_opts="-XX:ParallelGCThreads=10",
        mem_mb=4096
    wrapper:
        "v3.10.2/bio/snpeff/annotate"

rule running_multiqc:
    input:
        expand("annotated/{annotater}/{chr}.imputed.poly_subset.vcf.gz", chr=config['chromosomes'], annotater=['vep', 'snpeff']),
    output:
        "annotated/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "cd annotated && multiqc . --interactive"
     
