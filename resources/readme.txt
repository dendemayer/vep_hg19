# sample_names.tsv with:
zcat vcfs/chr1.imputed.poly_subset.vcf.gz | grep -v '##' | head -1 | cut -f 10- > resources/sample_names.tsv
