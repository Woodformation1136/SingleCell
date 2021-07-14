library(magrittr)
library(Matrix)
library(pheatmap)
library(RColorBrewer)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210525_Correlation_between_SCseq')


cal_cluster_mean_exp = function(UMI_matrix,cluster_df,species_name='Ptr',ncluster=10){
   foo = sapply(1:ncluster,function(cl){
      selected_barcode = with(cluster_df,Barcode[Cluster==cl])
      out = UMI_matrix %>% extract(,selected_barcode,drop=F) %>% apply(1,mean)
      return(out)
   }) %>% set_colnames(paste0(species_name,'_CellCluster_',1:ncluster))
   return(foo)
}

#Input data
##Input cluster data
Ptr_SCseq_cell_cluster_10 = readRDS('RDS_Ptr_SCseq_cell_cluster_10.rds')
Lch_SCseq_cell_cluster_10 = readRDS('RDS_Lch_SCseq_cell_cluster_10.rds')
Egr_SCseq_cell_cluster_10 = read.csv('Egr_clusters.csv')
Tar_SCseq_cell_cluster_18 = read.csv('Tar_clusters_18.csv')
Tar_SCseq_cell_cluster_18$Barcode %<>% sub('Kun','Ta',.)
Ptr_SCseq_cell_cluster_10$Cluster %>% table
Egr_SCseq_cell_cluster_10$Cluster %>% table
Tar_SCseq_cell_cluster_18$Cluster %>% table
Lch_SCseq_cell_cluster_10$Cluster %>% table

##Input Ptr data
selected_barcodes = Ptr_SCseq_cell_cluster_10$Barcode
Ptr_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Ptr_PtrEgr_UMI_matrix = readRDS('ortho_Ptr_Egr/RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Ptr_PtrTar_UMI_matrix = readRDS('ortho_Ptr_Tar/RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Ptr_PtrLch_UMI_matrix = readRDS('ortho_Ptr_Lch/RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Ptr_norm_factor = 1000/apply(Ptr_UMI_matrix,2,sum)
Ptr_PtrEgr_norm_UMI_matrix = apply(Ptr_PtrEgr_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
Ptr_PtrTar_norm_UMI_matrix = apply(Ptr_PtrTar_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
Ptr_PtrLch_norm_UMI_matrix = apply(Ptr_PtrLch_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix

Ptr_PtrEgr_mean_UMI_matrix = cal_cluster_mean_exp(Ptr_PtrEgr_UMI_matrix,Ptr_SCseq_cell_cluster_10,'Ptr')
Ptr_PtrTar_mean_UMI_matrix = cal_cluster_mean_exp(Ptr_PtrTar_UMI_matrix,Ptr_SCseq_cell_cluster_10,'Ptr')
Ptr_PtrLch_mean_UMI_matrix = cal_cluster_mean_exp(Ptr_PtrLch_UMI_matrix,Ptr_SCseq_cell_cluster_10,'Ptr')
Ptr_PtrEgr_mean_norm_UMI_matrix = cal_cluster_mean_exp(Ptr_PtrEgr_norm_UMI_matrix,Ptr_SCseq_cell_cluster_10,'Ptr')
Ptr_PtrTar_mean_norm_UMI_matrix = cal_cluster_mean_exp(Ptr_PtrTar_norm_UMI_matrix,Ptr_SCseq_cell_cluster_10,'Ptr')
Ptr_PtrLch_mean_norm_UMI_matrix = cal_cluster_mean_exp(Ptr_PtrLch_norm_UMI_matrix,Ptr_SCseq_cell_cluster_10,'Ptr')

##Input Egr data
selected_barcodes = Egr_SCseq_cell_cluster_10$Barcode
Egr_UMI_matrix = readRDS('RDS_Egr_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Egr_PtrEgr_UMI_matrix = readRDS('ortho_Ptr_Egr/RDS_Egr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Egr_norm_factor = 1000/apply(Egr_UMI_matrix,2,sum)
Egr_PtrEgr_norm_UMI_matrix = apply(Egr_PtrEgr_UMI_matrix,1,multiply_by,Egr_norm_factor) %>% t %>% Matrix

Egr_PtrEgr_mean_UMI_matrix = cal_cluster_mean_exp(Egr_PtrEgr_UMI_matrix,Egr_SCseq_cell_cluster_10,'Egr')
Egr_PtrEgr_mean_norm_UMI_matrix = cal_cluster_mean_exp(Egr_PtrEgr_norm_UMI_matrix,Egr_SCseq_cell_cluster_10,'Egr')

##Input Tar data
selected_barcodes = Tar_SCseq_cell_cluster_18$Barcode
Tar_UMI_matrix = readRDS('RDS_Tar_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Tar_PtrTar_UMI_matrix = readRDS('ortho_Ptr_Tar/RDS_Tar_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Tar_norm_factor = 1000/apply(Tar_UMI_matrix,2,sum)
Tar_PtrTar_norm_UMI_matrix = apply(Tar_PtrTar_UMI_matrix,1,multiply_by,Tar_norm_factor) %>% t %>% Matrix

Tar_PtrTar_mean_UMI_matrix = cal_cluster_mean_exp(Tar_PtrTar_UMI_matrix,Tar_SCseq_cell_cluster_18,'Tar',
                                                  ncluster=18)
Tar_PtrTar_mean_norm_UMI_matrix = cal_cluster_mean_exp(Tar_PtrTar_norm_UMI_matrix,Tar_SCseq_cell_cluster_18,'Tar',
                                                       ncluster=18)

##Input Lch data
selected_barcodes = Lch_SCseq_cell_cluster_10$Barcode
Lch_UMI_matrix = readRDS('RDS_Lch_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Lch_PtrLch_UMI_matrix = readRDS('ortho_Ptr_Lch/RDS_Lch_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
Lch_norm_factor = 1000/apply(Lch_UMI_matrix,2,sum)
Lch_PtrLch_norm_UMI_matrix = apply(Lch_PtrLch_UMI_matrix,1,multiply_by,Lch_norm_factor) %>% t %>% Matrix

Lch_PtrLch_mean_UMI_matrix = cal_cluster_mean_exp(Lch_PtrLch_UMI_matrix,Lch_SCseq_cell_cluster_10,'Lch')
Lch_PtrLch_mean_norm_UMI_matrix = cal_cluster_mean_exp(Lch_PtrLch_norm_UMI_matrix,Lch_SCseq_cell_cluster_10,'Lch')



plot_correlation_matrix = function(correlation_matrix = foo,
                                   rank_from_Ptr = F,
                                   plot_range = c(0,0.6),
                                   plot_width = 16,
                                   plot_height = 16,
                                   plot_name = 'PtrEgr_cluster_correlation.png',
                                   plot_legend = T,
                                   plot_axes = T){
   if(rank_from_Ptr){
      correlation_matrix = t(apply(correlation_matrix,1,function(x)rank(-x)))
      COL = colorRampPalette(brewer.pal(n=9,name='RdBu')[c(1,2,4,6,8,9)])(100)
      FPT = '%.d'
      CEX = 16
   }else{
      COL = colorRampPalette(rev(brewer.pal(n=9,name='RdBu')[c(1,2,5,8,9)]))(100)
      FPT = '%.4f'
      CEX = 12
   }
   png(plot_name,
       width=plot_width,height=plot_height,
       units='cm',res=300)
   pheatmap(correlation_matrix,cluster_rows=F,cluster_cols=F,scale='none',
            breaks=seq(plot_range[1],plot_range[2],length.out=100),
            color=COL,border_color='white',
            display_numbers=T,fontsize_number=CEX,number_format=FPT,
            show_rownames=plot_axes,show_colnames=plot_axes,legend=plot_legend)
   dev.off()
}


stopifnot(all(rownames(Ptr_PtrEgr_mean_norm_UMI_matrix)==rownames(Egr_PtrEgr_mean_norm_UMI_matrix)))
foo = cor(Ptr_PtrEgr_mean_norm_UMI_matrix[,c(1:8)],
          Egr_PtrEgr_mean_norm_UMI_matrix[,c(4,3,6,5,2,1,7,8)],method='pearson')
plot_correlation_matrix(plot_name = 'PtrEgr_cluster_correlation_whole.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 16,
                        plot_height = 16,
                        plot_legend = T,
                        plot_axes = T)
plot_correlation_matrix(plot_name = 'PtrEgr_cluster_correlation_legend.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 16,
                        plot_height = 16,
                        plot_legend = T,
                        plot_axes = F)
plot_correlation_matrix(plot_name = 'PtrEgr_cluster_correlation_clear.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 16,
                        plot_height = 16,
                        plot_legend = F,
                        plot_axes = F)
plot_correlation_matrix(plot_name = 'PtrEgr_cluster_correlation_rank.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = T,
                        plot_range = c(1,8),
                        plot_width = 16,
                        plot_height = 16,
                        plot_legend = F,
                        plot_axes = F)


stopifnot(all(rownames(Ptr_PtrLch_mean_norm_UMI_matrix)==rownames(Lch_PtrLch_mean_norm_UMI_matrix)))
foo = cor(Ptr_PtrLch_mean_norm_UMI_matrix[,c(1:8)],
          Lch_PtrLch_mean_norm_UMI_matrix[,c(1,2,5,8,4,6,3,7,9)])
plot_correlation_matrix(plot_name = 'PtrLch_cluster_correlation_whole.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = T,
                        plot_axes = T)
plot_correlation_matrix(plot_name = 'PtrLch_cluster_correlation_legend.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = T,
                        plot_axes = F)
plot_correlation_matrix(plot_name = 'PtrLch_cluster_correlation_clear.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = F,
                        plot_axes = F)
plot_correlation_matrix(plot_name = 'PtrLch_cluster_correlation_rank.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = T,
                        plot_range = c(1,9),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = F,
                        plot_axes = F)


stopifnot(all(rownames(Ptr_PtrTar_mean_norm_UMI_matrix)==rownames(Tar_PtrTar_mean_norm_UMI_matrix)))
foo = cor(Ptr_PtrTar_mean_norm_UMI_matrix[,c(1:8)],
          Tar_PtrTar_mean_norm_UMI_matrix[,c(1,4,6,9,8,5,2,7,3)])
plot_correlation_matrix(plot_name = 'PtrTar_cluster_correlation_whole.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = T,
                        plot_axes = T)
plot_correlation_matrix(plot_name = 'PtrTar_cluster_correlation_legend.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = T,
                        plot_axes = F)
plot_correlation_matrix(plot_name = 'PtrTar_cluster_correlation_clear.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = F,
                        plot_range = c(0,0.6),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = F,
                        plot_axes = F)
plot_correlation_matrix(plot_name = 'PtrTar_cluster_correlation_rank.png',
                        correlation_matrix = foo,
                        rank_from_Ptr = T,
                        plot_range = c(1,9),
                        plot_width = 18,
                        plot_height = 16,
                        plot_legend = F,
                        plot_axes = F)




Ptr_norm_UMI_matrix = apply(Ptr_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
Ptr_mean_norm_UMI_matrix = cal_cluster_mean_exp(Ptr_norm_UMI_matrix,Ptr_SCseq_cell_cluster_10,'Ptr')
COL = colorRampPalette(rev(brewer.pal(n=9,name='RdBu')[c(1,2,5,8,9)]))(100)
foo = cor(Ptr_mean_norm_UMI_matrix[,1:8])

png('Ptr_self_1_whole.png',
    width=16,height=16,
    units='cm',res=300)
pheatmap(foo,cluster_rows=F,cluster_cols=F,scale='none',
         breaks=seq(0,1,length.out=100),
         color=COL,border_color='white',
         display_numbers=T,fontsize_number=12,number_format='%.4f',
         show_rownames=T,show_colnames=T,legend=T)
dev.off()

png('Ptr_self_1_clear.png',
    width=16,height=16,
    units='cm',res=300)
pheatmap(foo,cluster_rows=F,cluster_cols=F,scale='none',
         breaks=seq(0,1,length.out=100),
         color=COL,border_color='white',
         display_numbers=T,fontsize_number=12,number_format='%.4f',
         show_rownames=F,show_colnames=F,legend=F)
dev.off()

png('Ptr_self_2_whole.png',
    width=16,height=16,
    units='cm',res=300)
pheatmap(foo[,8:1],cluster_rows=F,cluster_cols=F,scale='none',
         breaks=seq(0,1,length.out=100),
         color=COL,border_color='white',
         display_numbers=T,fontsize_number=12,number_format='%.4f',
         show_rownames=T,show_colnames=T,legend=T)
dev.off()

png('Ptr_self_2_clear.png',
    width=16,height=16,
    units='cm',res=300)
pheatmap(foo[,8:1],cluster_rows=F,cluster_cols=F,scale='none',
         breaks=seq(0,1,length.out=100),
         color=COL,border_color='white',
         display_numbers=T,fontsize_number=12,number_format='%.4f',
         show_rownames=F,show_colnames=F,legend=F)
dev.off()
