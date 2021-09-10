rule DEG_Pairwise:
    input:
        file = 'result/Mapping/Merged_seurat.rds'
    params:
        threads = 12
    output:
        txt = 'result/Pairwise/Pink1_venn_diagramm.pdf'
    shell:
        'Rscript workflow/scripts/2.1.DEG_Pairwise.R {input.file} {params.threads} {output.txt}'
