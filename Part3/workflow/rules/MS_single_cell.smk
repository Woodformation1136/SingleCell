# NOTE: Multi-sample single cell analysis -------------------------------------------
# Integrate SCseq samples
def get_UMI_csv_for_integration(wildcards) -> list:
    species_set = set([
        config["REF_SPECIES"][s.split("_")[1]]
        for s in config["MULTISAMPLES"][wildcards.multisample]
    ])
    if len(species_set) == 1:
        out = expand(
            "results/Single_species_analysis/all_UMI_tables/geneUMI_{sample}.csv",
            sample = config["MULTISAMPLES"][wildcards.multisample]
        )
    else:
        out = expand(
            "results/Ortholog_analysis/all_UMI_tables/orthologUMI_{sample}.csv",
            sample = config["MULTISAMPLES"][wildcards.multisample]
        )
    return out

rule integrate_sample:
    # Integration order is the same as the sample order within multi_UMI_csv
    input:
        multi_UMI_csv = get_UMI_csv_for_integration
    output:
        integration_rds = "results/Multi_species_analysis/all_data_rds/integration_{multisample}_base.rds"
    params:
        integrating_order = lambda wildcards: config["SEURAT_INTEGRATING_ORDER"][wildcards.multisample],
        k_params = lambda wildcards: config["SEURAT_K_PARAM"][wildcards.multisample],
        n_feature_cutoff_on_sample = 2000,
        n_barcode_cutoff_on_sample = 100
    script:
        "../scripts/integrate_scseq_data_with_CCA.R"

rule integrate_sample_runUMAP:
    conda:
        "../envs/PY39.yaml"
    input:
        integration_rds = "results/Multi_species_analysis/all_data_rds/integration_{multisample}_base.rds"
    output:
        integrated_umap_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_onlyUMAP_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    script:
        "../scripts/integrate_scseq_data_runUMAP.R"

# Turn Seurat CCA results into csv:
## 1. plotting_csv: plotting information of each barcode
## (columns: Barcode, Sample, Species, UMAP.1, UMAP.2, Cluster)
rule turn_CCA_results_into_csv:
    input:
        integrated_rds = "results/Multi_species_analysis/all_data_rds/integration_{multisample}_base.rds",
        integrated_umap_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_onlyUMAP_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    output:
        plotting_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    params:
        sample = lambda wildcards: config["MULTISAMPLES"][wildcards.multisample],
        species = lambda wildcards: [
            config["SPECIES"][sample]
            for sample in config["MULTISAMPLES"][wildcards.multisample]
        ]
    script:
        "../scripts/turn_seuratCCA_results_into_csv.R"

# Merge single sample (SS) and multi-sample (MS) csv
rule merge_SS_and_MS_plotting_csv:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables/plotting_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv",
        SS_plotting_csv_list = lambda wildcards: expand(
            "results/Single_species_analysis/all_plotting_tables/plotting_{sample}.csv",
            sample = config["MULTISAMPLES"][wildcards.multisample]
        )
    output:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    script:
        "../scripts/merge_SS_and_MS_plotting_csv.R"

# Plot integrated data by sample
rule plot_integrated_data_by_sample:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    output:
        figure_folder = directory("results/Multi_species_analysis/UMAP_by_sample/integration_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}")
    script:
        "../scripts/plot_integrated_data_by_samples.R"

# Plot integrated data by species
rule plot_integrated_data_by_species:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv"
    output:
        figure_folder = directory("results/Multi_species_analysis/UMAP_by_species/integration_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}")
    script:
        "../scripts/plot_integrated_data_by_species.R"

# Plot integrated data with correlation between MQD (LCM-seq) and SC-seq data
rule plot_integrated_data_with_correlation:
    input:
        MS_plotting_csv = "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}.csv",
        SS_MQD_plotting_csv = "results/Single_species_analysis/all_plotting_tables_cor_MQD/plotting_{SC_sample}_cor_{MQD_group}.csv"
    output:
        figure_folder = directory("results/Multi_species_analysis/UMAP_with_correlation/integration_{multisample}_seed_{seed_use}_md_{min_dist}_nn_{n_neighbors}/{SC_sample}_cor_{MQD_group}")
    script:
        "../scripts/plot_integrated_data_with_correlation.R"
