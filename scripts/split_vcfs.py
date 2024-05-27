import pandas as pd
import sys
import gzip

if "snakemake" in dir():
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")

    print('# snakemake inputs:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

    print('# snakemake output:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

    print('# snakemake wildcards:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

    print('# snakemake params:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.params.items()]

    vcf = snakemake.input.vcf
    vcf_out = snakemake.output.vcf_out
    sample = snakemake.wildcards.sample
    chr_ = snakemake.wildcards.chr

else:
    ################################################################################
    #                            3 projects test input                             #
    ################################################################################
    # snakemake inputs:
    vcf = "/scr/dings/PEVO/22_DNA_gencode_ref_new2/anne_vep/vcfs/chr1.imputed.poly_subset.vcf.gz"
    vcf_out = "/scr/dings/PEVO/22_DNA_gencode_ref_new2/anne_vep/vcfs_sample_split/chr1_DE72GSAHHMD100874.imputed.poly_subset.vcf.gz"
    sample = "DE72GSAHHMD100874"
    chr_ = "chr1"

count = 0
lines = []
with gzip.open(vcf, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            lines.append(line)
            count = count + 1
            continue
        else:
            break

with gzip.open(vcf_out, 'wt') as w:
    w.writelines(lines)

DF_index = pd.read_table(vcf, skiprows=count, nrows=100).iloc[:, 0:9]
DF_content = pd.read_table(vcf, usecols=[sample],skiprows=count, nrows=100)
print(f'saving {sample} specific vcf in {vcf_out}:')
pd.concat([DF_index, DF_content], axis=1).to_csv(vcf_out, compression='gzip', sep='\t', mode='a', index=False)

