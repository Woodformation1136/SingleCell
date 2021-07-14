library(magrittr)
library(dplyr)
library(slingshot)
library(pheatmap)
library(forecast)
library(viridis)
library(scales)
library(RColorBrewer)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\Single_cell\\20210603_PtrTar_Pseudotime_analysis_on_Seurat_UMAP')

PtrTar_allmarkers = readRDS('RDS_PtrTar_allmarkers.rds')
PtrTar_Seurat_projection_UMAP = readRDS('RDS_PtrTar_Seurat_projection_UMAP.rds')
Ptr_SCseq_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')
Tar_SCseq_UMI_matrix = readRDS('RDS_Tar_SCseq_UMI_matrix.rds')
PtrTar_Ptr_SCseq_ortho_UMI_matrix = readRDS('RDS_PtrTar_Ptr_SCseq_ortho_UMI_matrix.rds')
PtrTar_Tar_SCseq_ortho_UMI_matrix = readRDS('RDS_PtrTar_Tar_SCseq_ortho_UMI_matrix.rds')

Seurat_filter = function(UMI_matrix){
    features_per_cell = apply(UMI_matrix,2,function(x)sum(x>0))
    temp = UMI_matrix[,features_per_cell>=200]
    cells_per_feature = apply(temp,1,function(x)sum(x>0))
    out = temp[cells_per_feature>=3,]
    return(out)
}
PtrTar_Ptr_SCseq_ortho_UMI_matrix %<>% Seurat_filter
PtrTar_Tar_SCseq_ortho_UMI_matrix %<>% Seurat_filter

Ptr_selected_barcodes = sub('Ptr_','',colnames(PtrTar_Ptr_SCseq_ortho_UMI_matrix))
Tar_selected_barcodes = sub('Tar_','',colnames(PtrTar_Tar_SCseq_ortho_UMI_matrix))

Ptr_scale_factor = apply(Ptr_SCseq_UMI_matrix,2,function(x)1000/sum(x))[Ptr_selected_barcodes]
Tar_scale_factor = apply(Tar_SCseq_UMI_matrix,2,function(x)1000/sum(x))[Tar_selected_barcodes]

PtrTar_Ptr_SCseq_norm_ortho_UMI_matrix = apply(PtrTar_Ptr_SCseq_ortho_UMI_matrix,1,
                                               multiply_by,Ptr_scale_factor) %>% t
PtrTar_Tar_SCseq_norm_ortho_UMI_matrix = apply(PtrTar_Tar_SCseq_ortho_UMI_matrix,1,
                                               multiply_by,Tar_scale_factor) %>% t
PtrTar_overlap_ortholog = intersect(rownames(PtrTar_Ptr_SCseq_norm_ortho_UMI_matrix),
                                    rownames(PtrTar_Tar_SCseq_norm_ortho_UMI_matrix))
PtrTar_SCseq_norm_ortho_UMI_matrix = cbind(PtrTar_Ptr_SCseq_norm_ortho_UMI_matrix[PtrTar_overlap_ortholog,],
                                           PtrTar_Tar_SCseq_norm_ortho_UMI_matrix[PtrTar_overlap_ortholog,])


#Identify the psuedotime of each cell on each lineage
PtrTar_lineage = new('SlingshotDataSet')
str(PtrTar_lineage)

##Prepare the reducedDim
rd = as.matrix(PtrTar_Seurat_projection_UMAP[,1:2])

##Prepare the clusterLabels
cluster_factor = factor(PtrTar_Seurat_projection_UMAP$Seurat_color)
levels(cluster_factor)
cl = sapply(as.numeric(cluster_factor),
            function(i){
                out = rep(0,10); out[i] = 1
                return(out)
            }) %>% t
rownames(cl) = rownames(rd)
colnames(cl) = 1:10

##Prepare the lineages
lin = list(Lineage1 = c('9','8','7'),
           Lineage2 = c('4','2','3','1'),
           Lineage3 = c('4','2','3','6'))

##Prepare the adjacency
adc = matrix(0,10,10)
rownames(adc) = 1:10
colnames(adc) = 1:10
adc[9,8] = adc[8,9] = 1
adc[8,7] = adc[7,8] = 1
adc[4,2] = adc[2,4] = 1
adc[2,3] = adc[3,2] = 1
adc[3,1] = adc[1,3] = 1
adc[3,6] = adc[6,3] = 1

##Fill in the contents
PtrTar_lineage@reducedDim = rd
PtrTar_lineage@clusterLabels = cl
PtrTar_lineage@lineages = lin
PtrTar_lineage@adjacency = adc

##Plot the lineages
plot(PtrTar_Seurat_projection_UMAP$UMAP_1,
     PtrTar_Seurat_projection_UMAP$UMAP_2,
     pch=20,cex=0.3,
     col=PtrTar_Seurat_projection_UMAP$Seurat_color,
     xlab='UMAP_1',ylab='UMAP_2',main='PtrTar')
lines(PtrTar_lineage, lwd=3, col='black')

##Get lineage curves
PtrTar_lineage = getCurves(PtrTar_lineage)
plot(PtrTar_Seurat_projection_UMAP$UMAP_1,
     PtrTar_Seurat_projection_UMAP$UMAP_2,
     pch=20,cex=0.3,
     col=PtrTar_Seurat_projection_UMAP$Seurat_color,
     xlab='UMAP_1',ylab='UMAP_2',main='PtrTar')
lines(PtrTar_lineage, lwd=3, col='black')

##Get psuedotime matrix
PtrTar_psuedotime_matrix = slingPseudotime(PtrTar_lineage)
PtrTar_lineage_list = sapply(1:3,
                             function(lid){
                                 lineage_barcodes = PtrTar_psuedotime_matrix[,lid] %>% 
                                     extract(!is.na(.)) %>% sort %>% names
                             })
names(PtrTar_lineage_list) = c('Ray','Vessel','Fiber')
PtrTar_total_lineage = c(PtrTar_lineage_list$Ray,
                         PtrTar_lineage_list$Vessel,
                         PtrTar_lineage_list$Fiber)

##Plot each lineage
png('PtrTar_Psuedotime_lineages.png',
    pointsize=10,width=60,height=15,units='cm',res=300)
par(mfcol=c(1,3),mar=c(0,0,0,0))
COL = rev(inferno(150)[26:125])
for(lineage_id in 1:3){
    plot_psuedotime = na.omit(PtrTar_psuedotime_matrix[,lineage_id])
    col_ind = round(plot_psuedotime/max(plot_psuedotime)*99)+1
    plot_order = match(names(col_ind),
                       PtrTar_Seurat_projection_UMAP$SpBarcode)
    plot(PtrTar_Seurat_projection_UMAP$UMAP_1,
         PtrTar_Seurat_projection_UMAP$UMAP_2,
         pch=20,col='#EEEEEE',cex=2,axes=F,
         xlab='',ylab='')
    legend('topleft',c('Ray','Vessel','Fiber')[lineage_id],bty='n',cex=3)
    points(PtrTar_Seurat_projection_UMAP$UMAP_1[plot_order],
           PtrTar_Seurat_projection_UMAP$UMAP_2[plot_order],
           col=COL[col_ind],pch=20,cex=3)
}
par(oriPar)
dev.off()

png('PtrTar_Psuedotime_lineages_legend.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(seq(0,1,length.out=length(COL)),
     rep(0,length(COL)),
     col=COL,pch='|',cex=7)
dev.off()


heatmap_based_on_lineage = function(plot_UMI_matrix = PtrTar_SCseq_norm_ortho_UMI_matrix,
                                    plot_species = 'Both',
                                    plot_lineage = 'Total',
                                    plot_features_list = Total_DE_ortholog,
                                    plot_name = 'Rplot.png',
                                    ma_order = 50,
                                    output_figure = T){
    if(plot_lineage=='Total'){
        plot_barcodes_list = PtrTar_lineage_list
        if(plot_species!='Both')
            plot_barcodes_list %<>% sapply(function(x)x[grepl(plot_species,x)])
        lineage_length = sapply(plot_barcodes_list,length)
        plot_barcodes = unlist(plot_barcodes_list)
    }else{
        plot_barcodes = PtrTar_lineage_list[[plot_lineage]]
        if(plot_species!='Both')
            plot_barcodes %<>% extract(grepl(plot_species,.))
        lineage_length = length(plot_barcodes)
    }
    
    plot_features = unlist(plot_features_list)
    if(length(lineage_length)==1){
        plot_expression_matrix = plot_UMI_matrix[plot_features,plot_barcodes] %>%
            apply(1,ma,ma_order) %>% add(1) %>% log2 %>% t
    }else{
        lineage_length_matrix = matrix(c(head(c(0,cumsum(lineage_length))+1,-1),
                                         tail(c(0,cumsum(lineage_length)),-1)),
                                       nrow=2,byrow=T)
        plot_expression_matrix = plot_UMI_matrix[plot_features,plot_barcodes] %>%
            apply(1,function(x)
                unlist(apply(lineage_length_matrix,2,
                             function(y) ma(x[y[1]:y[2]],ma_order)))) %>% 
            add(1) %>% log2 %>% t
    }
    
    temp_barcodes = paste0('Cell_',seq(ncol(plot_expression_matrix)))
    colnames(plot_expression_matrix) = temp_barcodes
    temp_features = paste0('Feature_',seq(nrow(plot_expression_matrix)))
    rownames(plot_expression_matrix) = temp_features
    
    # print_as_Cluster_05d = function(cluster_id){
    #     out = cluster_id %>% sub('Cluster_','',.) %>%
    #         as.numeric %>% sprintf('%05d',.) %>% paste0('Cluster_',.)
    #     return(out)
    # }
    # rownames(plot_expression_matrix) %<>% print_as_Cluster_05d
    # plot_unique_features = unique(plot_features)
    # plot_mean_UMI = apply(plot_UMI_matrix[plot_unique_features,],1,mean)
    # plot_max_UMI = apply(plot_UMI_matrix[plot_unique_features,],1,max)
    # anno_row = data.frame(MeanUMI=plot_mean_UMI,
    #                       MaxUMI=plot_max_UMI,
    #                       row.names=temp_features)
    
    feature_group = rep(head(LETTERS,length(plot_features_list)),
                        times=sapply(plot_features_list,length))
    plot_mean_UMI = apply(plot_UMI_matrix[plot_features,],1,mean)
    plot_max_UMI = apply(plot_UMI_matrix[plot_features,],1,max)
    anno_row = data.frame(MeanUMI=plot_mean_UMI,
                          MaxUMI=plot_max_UMI,
                          FeatureGroup = feature_group,
                          row.names=temp_features)
    
    if(plot_species=='Both'){
        anno_col = data.frame(Cluster=PtrTar_Seurat_projection_UMAP[plot_barcodes,'Seurat_color'],
                              Species=substr(plot_barcodes,1,3),
                              row.names=temp_barcodes)
        ann_colors = list(Cluster=set_names(unique(anno_col[,1]),
                                            unique(anno_col[,1])),
                          # Species=c('Ptr'=hue_pal()(2)[1],
                          #           'Tar'=hue_pal()(2)[2]),
                          Species=c('Ptr'='black',
                                    'Tar'='#C59739'),
                          FeatureGroup=set_names(brewer.pal(length(plot_features_list),'Dark2'),
                                                 head(LETTERS,length(plot_features_list))))
    }else{
        anno_col = data.frame(Cluster=PtrTar_Seurat_projection_UMAP[plot_barcodes,'Seurat_color'],
                              row.names=temp_barcodes)
        ann_colors = list(Cluster=set_names(unique(anno_col[,1]),
                                            unique(anno_col[,1])),
                          FeatureGroup=set_names(brewer.pal(length(plot_features_list),'Dark2'),
                                                 head(LETTERS,length(plot_features_list))))
    }
    
    
    scale_plot_expression_matrix = plot_expression_matrix %>% 
        # apply(1,scale,center=T,scale=T) %>%
        # apply(1,function(x)100*x/max(x,na.rm=T)) %>%
        apply(1,function(x)100*punif(x,min(x,na.rm=T),max(x,na.rm=T))) %>%
        t %>% set_colnames(temp_barcodes)
    
    
    scale_plot_expression_matrix %<>% extract(,!is.na(.[1,]))
    lineage_length %<>% subtract(ma_order - (ma_order%%2))
    
    out_heatmap = pheatmap(scale_plot_expression_matrix,
                           clustering_method = 'ward.D2',
                           clustering_distance_rows = 'manhattan',
                           cluster_rows = F,
                           cluster_cols = F,
                           color = viridis(100),
                           breaks = seq(0,max(scale_plot_expression_matrix,na.rm=T),
                                        length.out=100),
                           border_color = viridis(100)[1],
                           gaps_row = cumsum(sapply(plot_features_list,length)),
                           gaps_col = cumsum(lineage_length),
                           show_rownames = F,
                           show_colnames = F,
                           fontsize_row = 8,
                           annotation_row = anno_row,
                           annotation_col = anno_col,
                           annotation_colors = ann_colors,
                           silent = output_figure,
                           na_col = 'gray80')
    
    if(output_figure){
        png(plot_name, width=40, height=20, units='cm',res=600)
        grid::grid.newpage()
        grid::grid.draw(out_heatmap$gtable)
        dev.off()
    }
}



Ray_lineage_color = c('#FF7F0F','#F8E71C','#E377C2')
Vessel_lineage_color = c('#9467BD','#2AA02A','#8C564B','#1F77B4')
Fiber_lineage_color = c('#9467BD','#2AA02A','#8C564B','#D62728')
Total_lineage_color = c(Ray_lineage_color,Vessel_lineage_color,Fiber_lineage_color[4])
Total_DE_ortholog = lapply(Total_lineage_color,
                           function(COL){
                               out = PtrTar_allmarkers %>%
                                   filter(p_val_adj<0.05,cluster==COL) %>%
                                   arrange(-avg_log2FC)
                               out$gene %<>% sub('-','_',.)
                               return(out)
                           })
for(i in seq_along(Total_DE_ortholog)){
    write.csv(Total_DE_ortholog[[i]],
              file=paste0(LETTERS[i],'.csv'),row.names=F)
}

heatmap_based_on_lineage(plot_UMI_matrix = PtrTar_SCseq_norm_ortho_UMI_matrix,
                         plot_species = 'Both',
                         plot_lineage = 'Total',
                         plot_features_list = Total_DE_ortholog,
                         plot_name = 'PtrTar_Both_DE_ortholog_heatmap.png',
                         ma_order = 21,
                         output_figure = T)

heatmap_based_on_lineage(plot_UMI_matrix = PtrTar_SCseq_norm_ortho_UMI_matrix,
                         plot_species = 'Ptr',
                         plot_lineage = 'Total',
                         plot_features_list = Total_DE_ortholog,
                         plot_name = 'PtrTar_Ptr_DE_ortholog_heatmap.png',
                         ma_order = 21,
                         output_figure = T)

heatmap_based_on_lineage(plot_UMI_matrix = PtrTar_SCseq_norm_ortho_UMI_matrix,
                         plot_species = 'Tar',
                         plot_lineage = 'Total',
                         plot_features_list = Total_DE_ortholog,
                         plot_name = 'PtrTar_Tar_DE_ortholog_heatmap.png',
                         ma_order = 21,
                         output_figure = T)

