# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

def get_dem(wildcards):
    # no trimming, use raw reads
    return samples.loc[(wildcards.sample), ["file"]]



# file=expand('DATA/{file}',file=FILES)
rule Preprocess:
    input:
        file=get_dem
    params:
        sample='{sample}'
    output:
        txt = 'result/Preprocess/{sample}/{sample}_{condition}_Cells_Removed.rds'
    shell:
        'Rscript workflow/scripts/0.Preprocessing.R {input.file} {params.sample} {output.txt}'

# txt=expand("result/{sample}_Cells_Removed.tsv",sample=samples["sample"])
