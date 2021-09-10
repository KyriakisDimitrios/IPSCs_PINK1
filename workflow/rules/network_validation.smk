rule Network_Validation:
    input:
        file = 'result/Mapping/Merged_seurat.rds',
        #groupA = "result/DEG/DF_Group_A_Conserved_all.txt",
        groupB = "result/DEG/DF_Group_B_Conserved_all.txt",
        groupC = "result/DEG/DF_Group_C_Conserved_all.txt",
        groupD = "result/DEG/DF_Group_D.txt"
    params:
        threads = 4
    output:
        "result/Network/Num_Nodes_Interactions.pdf"
    shell:
        "Rscript workflow/scripts/4.Network_Validation.R {input.file} {input.groupB} {input.groupC} {input.groupD} {output}"
