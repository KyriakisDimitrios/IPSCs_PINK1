rule DEG_Group_BC:
    input:
        file = 'result/Mapping/Merged_seurat.rds'
    params:
        threads = 4
    output:
        groupB = "result/DEG/DF_Group_B_Conserved_all.txt",
        groupC = "result/DEG/DF_Group_C_Conserved_all.txt"
    shell:
        "Rscript workflow/scripts/2.2.DEG_Group_BC.R {input.file} {params.threads} {output.groupB} {output.groupC}"
