def get_rds(wildcards):
    files_list=expand('result/Preprocess/{sample}/{sample}_{condition}_Cells_Removed.rds',zip,sample=samples['sample'],condition=samples['Condition'])
    return files_list



rule Mapping:
    input:
        file = get_rds,
    params:
        threads = 6
    output:
        txt = 'result/Mapping/Merged_seurat.rds'
    shell:
        'Rscript workflow/scripts/1.Mapping.R {input.file} {params.threads} {output.txt}'

# txt=expand("result/{sample}_Cells_Removed.tsv",sample=samples["sample"])
