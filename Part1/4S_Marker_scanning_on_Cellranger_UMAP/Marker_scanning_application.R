library(magrittr)
library(Matrix)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210521_Marker_scanning')

load('Plot_specific_feature.RData')

#Usage:
# plot_specific_feature(Species = c('Ptr','Egr','Tar','Lch'),
#                       species_group = c('Ptr','PtrEgr','PtrTar','PtrLch','4S'),
#                       Norm = c(T,F),
#                       point_size = 1,
#                       featureID = c('Potri.004G059600.v4.1','Cluster_0'),
#                       summary_title = 'Summary UMI',
#                       summary_method = c('sum','max'),
#                       output_filename = 'Ptr_CellRanger_UMAP.png',
#                       color_scale = c(20,50),
#                       output_for_Jimmy = c(F,T))
# plot_specific_feature_4in1(Norm = c(T,F),
#                            point_size = 1,
#                            featureID = 'Cluster_0',
#                            summary_title = 'Summary UMI',
#                            summary_method = c('sum','max'),
#                            output_filename = '4S_CellRanger_UMAP.png',
#                            color_scale = list(Ptr=c(20,50),Egr=c(20,50),
#                                               Tar=c(20,50),Lch=c(20,50)),
#                            output_for_Jimmy = F)
#Notice:
#1. If length(featureID)==1, ignore summary_title and summary_method
#2. If summary_method!='none', summary_title can be assigned as better plot title
#3. Set color_scale[1] < color_scale[2], and both between 0 to 100
#4. If the result is HOW-KAN--, set output_for_Jimmy=T and send it to Jimmy
#
#Examples:
#1. Single species plot
#(1) Output UMAP of Ptr with normlized UMI of geneID Potri.004G059600.v4.1
#    using color turning percentile as 20 and 50
plot_specific_feature(Species = 'Ptr',
                      species_group = 'Ptr',
                      Norm = T,
                      featureID = 'Potri.004G059600.v4.1',
                      output_filename = 'Ptr_CellRanger_UMAP.png',
                      color_scale = c(20,50))
#(2) Output UMAP of Egr with normlized UMI of geneID Eucgr.A00001.v2.0
#    using color turning percentile as 20 and 50
plot_specific_feature(Species = 'Egr',
                      species_group = 'Egr',
                      Norm = T,
                      featureID = 'Eucgr.A00001.v2.0',
                      output_filename = 'Ptr_CellRanger_UMAP.png',
                      color_scale = c(20,50))
#(3) Output UMAP of Ptr with non-normalized UMI of geneID Potri.004G059600.v4.1
#    using color turning percentile as 40 and 80
plot_specific_feature(Species = 'Ptr',
                      species_group = 'Ptr',
                      Norm = F,
                      featureID = 'Potri.004G059600.v4.1',
                      output_filename = 'Ptr_CellRanger_UMAP_norm.png',
                      color_scale = c(40,80))
#(4) Output UMAP of Ptr with total normlized UMI of geneIDs
#    Potri.004G059600.v4.1 and Potri.006G181900.v4.1
#    using color turning percentile as 20 and 50
plot_specific_feature(Species = 'Ptr',
                      species_group = 'Ptr',
                      Norm = T,
                      featureID = c('Potri.004G059600.v4.1','Potri.006G181900.v4.1'),
                      summary_title = 'Two for Example',
                      summary_method = 'sum',
                      output_filename = 'Ptr_CellRanger_UMAP_total.png',
                      color_scale = c(20,50))
#(5) Output UMAP of Ptr with normalized UMI of clusterID Cluster_0
#    comparing the ortholog between Ptr and Egr
#    using color turning percentile as 40 and 80
plot_specific_feature(Species = 'Ptr',
                      species_group = 'PtrEgr',
                      Norm = T,
                      featureID = 'Cluster_0',
                      output_filename = 'Ptr_PtrEgr_CellRanger_UMAP.png',
                      color_scale = c(40,80))
#(6) Output UMAP of Tar with non-normalized UMI of clusterID Cluster_0
#    comparing the ortholog among 4 species
#    using color turning percentile as 20 and 50
plot_specific_feature(Species = 'Tar',
                      species_group = '4S',
                      Norm = F,
                      featureID = 'Cluster_0',
                      output_filename = 'Tar_PtrTar_CellRanger_UMAP.png',
                      color_scale = c(20,50))
#2. Four species in one plot
#(1) Output UMAP of 4 species with normalized UMI of clusterID Cluster_0
#    using color turning percentile as 20 and 50 for each species
plot_specific_feature_4in1(Norm = T,
                           featureID = 'Cluster_0',
                           output_filename = '4S_CellRanger_UMAP.png',
                           color_scale = list(Ptr=c(20,50),Egr=c(20,50),
                                              Tar=c(20,50),Lch=c(20,50)),
                           output_for_Jimmy = F)
#(2) Output UMAP of 4 species with normalized UMI of clusterID Cluster_0
#    using color turning percentile as 5 and 35 for Ptr,
#    10 and 40 for Egr, 15 and 45 for Tar, 20 and 50 for Lch,
plot_specific_feature_4in1(Norm = T,
                           featureID = 'Cluster_0',
                           output_filename = '4S_CellRanger_UMAP.png',
                           color_scale = list(Ptr=c(5,35),Egr=c(10,40),
                                              Tar=c(15,45),Lch=c(20,50)),
                           output_for_Jimmy = F)
#(3) Output UMAP of 4 species with maximum non-normlized UMI between
#    clusterIDs Cluster_0 and Cluster_1
#    using color turning percentile as 20 and 50 for each species
plot_specific_feature_4in1(Norm = F,
                           featureID = c('Cluster_0','Cluster_1'),
                           summary_title = 'Summary UMI',
                           summary_method = 'max',
                           output_filename = '4S_CellRanger_UMAP.png',
                           color_scale = list(Ptr=c(20,50),Egr=c(20,50),
                                              Tar=c(20,50),Lch=c(20,50)),
                           output_for_Jimmy = F)
