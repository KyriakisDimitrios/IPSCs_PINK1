rule DEG_Group_D:
    input:
        file = 'result/Mapping/Merged_seurat.rds'
    params:
        threads = 4
    output:
        groupD = "result/DEG/DF_Group_D.txt"
    shell:
        "Rscript workflow/scripts/2.3.DEG_Group_D.R {input.file} {params.threads} {output.groupD}"
