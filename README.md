# VEP
 
## this Snakemake workflow initializes all needed environments and reference files (hg19) to annotate vcfs with Ensembl Variant Effect Predictor (VEP) and SNPeff:  

### example multiqc in annotated/multiqc_report.html

# Start the Workflow:  
within the repo issue:  
```bash
snakemake -p --cores 40 --use-conda                                                                                                               ±[●●][master]
```

the vcfs have to be placed within the vcfs dir like:  
```
vcfs/chr1.imputed.poly_subset.vcf.gz
```
depending on what vcfs are available change that within the config.yaml and add the chromosomes:
within config.yaml:  
```bash
chromosomes:
  - chr1
  - chr2
```
of course, the chr2 vcf must be present then too like so:  
```
vcfs/chr1.imputed.poly_subset.vcf.gz
vcfs/chr2.imputed.poly_subset.vcf.gz
```

multiqc aggregates everything applied within the vcfs dir and the config.yaml


VEP:  
    - https://www.ensembl.org/info/docs/tools/vep/index.html

reference fa used:  
    - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

gff3 used:
    - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz

snakemake wrapper used:
    - https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/vep/annotate.html
    - 

