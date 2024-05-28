configfile: 'config.yaml'

chr_sample_combinations = [i + '_' + j for i in  config['chromosomes'] for j in config['samples']]
rule all:
    input:
        "annotated/vcfs/multiqc_report.html",
        # requesting this multiqc triggers about 1400 jobs:
        #"annotated/vcfs_sample_split/multiqc_report.html"
        # those are the peptid list filtered annotated vcfs:
        expand("annotated/vcfs/{annotater}/{chr}.imputed.poly_subset_peptid_filtered.vcf.gz", annotater=['vep', 'snpeff'], chr=config['chromosomes']),
        # additional filters can be applied to the annotated vcfs:
        expand("annotated/vcfs/{annotater}/{chr}.imputed.poly_subset_peptid_filtered_protein_coding.vcf.gz", annotater=['vep', 'snpeff'], chr=config['chromosomes']),

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

rule limit_gff_to_gene:
    input:
        "ref_seq/{annot_file}_sorted.gff3",
    output:
        "ref_seq/{annot_file}_sorted_gene.gff3",
    shell:
        "grep -w gene {input} > {output}"

rule limit_gff_to_peptid:
    input:
        gff = "ref_seq/{annot_file}_sorted_gene.gff3",
        peptid = "resources/Peptide-receptor-gut_Tabelle1.tsv"
    output:
        "ref_seq/{annot_file}_sorted_gene_peptid.gff3"
    shell:
        "grep -f {input.peptid} {input.gff} > {output}"

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
        calls="{vcf_type}/{chr}.imputed.poly_subset.vcf.gz",  # .vcf, .vcf.gz or .bcf
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
        calls="annotated/{vcf_type}/vep/{chr}.imputed.poly_subset.vcf.gz",  # .vcf, .vcf.gz or .bcf stats="vep/variants.html",
        #--stats_file [SAMPLE_NAME].vep.html (without the vep or summary suffix, MultiQC will ignore the HTML files)
        stats='annotated/{vcf_type}/vep/{chr}.imputed.poly_subset.vep.html',
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra="--everything",  # optional: extra arguments
    benchmark: 
        "benchmarks/vep/{vcf_type}/{chr}_benchmark.log"
    log:
        "logs/vep/{vcf_type}_{chr}.log",
    threads: 3
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
        calls="{vcf_type}/{chr}.imputed.poly_subset.vcf.gz",  # .vcf, .vcf.gz or .bcf
        db="resources/snpeff/GRCh38.p14" # path to reference db downloaded with the snpeff download wrapper
    output:
        stats="annotated/{vcf_type}/snpeff/{chr}.imputed.poly_subset.html",  # summary statistics (in HTML), optional
        calls="annotated/{vcf_type}/snpeff/{chr}.imputed.poly_subset.vcf.gz",   # annotated calls (vcf, bcf, or vcf.gz)
        csvstats="annotated/{vcf_type}/snpeff/{chr}.imputed.poly_subset.csv", # summary statistics in CSV, optional
    log:
        "logs/snpeff/{vcf_type}_{chr}.log"
    resources:
        java_opts="-XX:ParallelGCThreads=10",
        mem_mb=4096
    threads:
        3
    wrapper:
        "v3.10.2/bio/snpeff/annotate"

rule running_multiqc:
    input:
        expand("annotated/vcfs/{annotater}/{chr}.imputed.poly_subset.vcf.gz", chr=config['chromosomes'], annotater=['vep', 'snpeff']),
    output:
        "annotated/vcfs/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "cd annotated/vcfs && multiqc . --interactive -f"

rule running_multiqc_split:
    input:
        expand("annotated/vcfs_sample_split/{annotater}/{chr}.imputed.poly_subset.vcf.gz", chr=chr_sample_combinations, annotater=['vep', 'snpeff']),
    output:
        "annotated/vcfs_sample_split/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "cd annotated/vcfs_sample_split && multiqc . --interactive -f"
     
rule split_vcfs:
    input:
        vcf = "vcfs/{chr}.imputed.poly_subset.vcf.gz",
    output:
        vcf_out = "vcfs_sample_split/{chr}_{sample}.imputed.poly_subset.vcf.gz"
    conda:
        "envs/pandas.yaml"
    log:
        "logs/split_vcfs/{chr}_{sample}.log"
    script:
        "scripts/split_vcfs.py"

rule filter_peptid_table:
    """
    filter the annotated vcfs on the given peptid table, works for both, snpeff and vep
    additionally filter on protein_coding flag within annotated vcfs:
    """
    input:
        vcf = "annotated/{vcf_type}/{annotater}/{chr}.imputed.poly_subset.vcf.gz",  # .vcf, .vcf.gz or .bcf stats="vep/variants.html",
        peptid_table = "resources/Peptide-receptor-gut_Tabelle1.tsv"
    output:
        "annotated/{vcf_type}/{annotater}/{chr}.imputed.poly_subset_peptid_filtered.vcf.gz",  # .vcf, .vcf.gz or .bcf stats="vep/variants.html",
    conda:
        "envs/bcftools.yaml"
    params:
        temp = "annotated/{vcf_type}/{annotater}/{chr}.temp"
    shell:
        """
        fgrep -iwf {input.peptid_table} <(bcftools view -H {input.vcf}) > {params.temp} || true
        cat <(bcftools view -h {input.vcf}) {params.temp} | gzip > {output}
        rm -r {params.temp}
        """

rule filter_annotation_flag:
    input:
        "annotated/{vcf_type}/{annotater}/{chr}.imputed.poly_subset_peptid_filtered.vcf.gz",
    output:
        "annotated/{vcf_type}/{annotater}/{chr}.imputed.poly_subset_peptid_filtered_{flag}.vcf.gz",
    conda:
        "envs/bcftools.yaml"
    params:
        temp = "annotated/{vcf_type}/{annotater}/{chr}_{flag}.temp"
    shell:
        """
        fgrep -iw {wildcards.flag} <(bcftools view -H {input}) > {params.temp} || true
        cat <(bcftools view -h {input}) {params.temp} | gzip > {output}
        rm -r {params.temp}
        """
