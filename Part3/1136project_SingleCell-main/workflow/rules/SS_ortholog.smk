# NOTE: Orthlog analysis -------------------------------------------------------------
# Create protein id to gene id table (Ptr, Egr, Tar, Lch) (customized process in R code)
rule get_protein2gene_table:
    params:
        Ptr_protein_fa = "/home/woodydrylab/HDD/GenomicsData/Populus_trichocarpa/Ptrichocarpa_533_v4.1/Ptrichocarpa_533_v4.1.protein_primaryTranscriptOnly.fa",
        Egr_gff = "/home/woodydrylab/HDD/GenomicsData/Eucalyptus_grandis/Egrandis_297_v2.0/Egrandis_297_v2.0.gene.gff3",
        Tar_gtf = "/home/woodydrylab/HDD/GenomicsData/Trochodendron_aralioides/Taralioides_20200702/Trochodendron_aralioides_chromosomes_pasa2.longest.filter.gtf",
        Lch_gff = "/home/woodydrylab/HDD/GenomicsData/Liriodendron_chinense/Liriodendron_chinense_20191212/Liriodendron_chinense_gff"
    output:
        "results/protein2gene.csv"
    script:
        "../scripts/get_protein2gene.R"

# Preprocess ortholog data
## 1. melt ortholog table from wide to long (for orthoMCL)
## 2. turn protein id into gene id
rule preprocess_ortho_data:
    input:
        orthoMCL_txt = lambda wildcards: config["RAWDATA_PATH"] + "/20220000_Single_cell_rawdata/20210130_Geneclustering_OrthoMCL/all.group",
        protein2gene_csv = "results/protein2gene.csv"
    output:
        ortholog_long_csv = "results/Ortholog_analysis/all_group_long.csv",
        ortholog_long_convertedID_csv = "results/Ortholog_analysis/all_group_long_convertedID.csv"
    script:
        "../scripts/preprocess_ortholog_data.R"

# Filter long ortholog table by getting unique transcript for each gene
rule get_primary_ortho_data:
    input:
        "results/Ortholog_analysis/all_group_long.csv"
    output:
        "results/Ortholog_analysis/primary_group_long.csv"
    params:
        rawdata_path = config["RAWDATA_PATH"]
    script:
        "../scripts/get_primary_ortholog.R"

# Calculate ortholog UMI count
rule compute_ortholog_expression:
    input:
        ortholog_long_convertedID_csv = "results/Ortholog_analysis/all_group_long_convertedID.csv",
        geneUMI_csv = "results/Single_species_analysis/all_UMI_tables/geneUMI_{sample}.csv"
    output:
        orthologUMI_csv = "results/Ortholog_analysis/all_UMI_tables/orthologUMI_{sample}.csv"
    script:
        "../scripts/calculate_ortholog_UMI_counts.R"