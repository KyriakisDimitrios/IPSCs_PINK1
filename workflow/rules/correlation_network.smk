rule correlation_network:
    input:
        file = 'result/Mapping/Merged_seurat.rds',
        csv1 = 'DATA/Networks/NODES_and_pathways_11.6.20_9.csv',
        csv2 = 'DATA/Networks/RUN_12_EDGES_ST_GM_REMERGE_used_11_6_20.csv',
        csv3 = 'DATA/Networks/NODES_May_29_Manual_Mito_Ubiq.csv',
        csv4 = 'DATA/Networks/EDGES_part1_manual_May_29.csv'
    output:
        'Correlation_Networks.rds'
    shell:
        'Rscript workflow/scripts/4.Correlation_Network.R {input.file} {input.csv1} {input.csv2} {input.csv3} {input.csv4} {output}'


