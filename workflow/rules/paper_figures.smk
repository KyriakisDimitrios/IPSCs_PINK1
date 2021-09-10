rule Paper_Figures:
    input:
        file = 'result/Mapping/Merged_seurat.rds'
    params:
        threads = 2
    output:
        txt = 'result/Paper_Figures/Figure4.pdf'
    shell:
        'Rscript workflow/scripts/3.Paper_Plots.R {input.file} {params.threads} {output.txt}'
