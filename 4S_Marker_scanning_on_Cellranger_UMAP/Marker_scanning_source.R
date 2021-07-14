library(magrittr)
library(Matrix)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210521_Marker_scanning')

# #Input data
# ##Input Ptr data
# Ptr_projection_umap = readRDS('RDS_Ptr_projection_umap.rds')
# selected_barcodes = Ptr_projection_umap$Barcode
# 
# Ptr_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Ptr_4S_UMI_matrix = readRDS('ortho_4S/RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Ptr_PtrEgr_UMI_matrix = readRDS('ortho_Ptr_Egr/RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Ptr_PtrTar_UMI_matrix = readRDS('ortho_Ptr_Tar/RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Ptr_PtrLch_UMI_matrix = readRDS('ortho_Ptr_Lch/RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Ptr_PtrUnion3S_UMI_matrix = Ptr_PtrEgr_UMI_matrix %>%
#     rbind(.,Ptr_PtrTar_UMI_matrix[setdiff(rownames(Ptr_PtrTar_UMI_matrix),
#                                           rownames(.)),]) %>%
#     rbind(.,Ptr_PtrLch_UMI_matrix[setdiff(rownames(Ptr_PtrLch_UMI_matrix),
#                                           rownames(.)),])
# rm(list=c('Ptr_PtrEgr_UMI_matrix','Ptr_PtrTar_UMI_matrix','Ptr_PtrLch_UMI_matrix')); gc()
# 
# Ptr_norm_factor = 1000/apply(Ptr_UMI_matrix,2,sum)
# Ptr_norm_UMI_matrix = apply(Ptr_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
# Ptr_4S_norm_UMI_matrix = apply(Ptr_4S_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
# Ptr_PtrUnion3S_norm_UMI_matrix = apply(Ptr_PtrUnion3S_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
# # Ptr_PtrEgr_norm_UMI_matrix = apply(Ptr_PtrEgr_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
# # Ptr_PtrTar_norm_UMI_matrix = apply(Ptr_PtrTar_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
# # Ptr_PtrLch_norm_UMI_matrix = apply(Ptr_PtrLch_UMI_matrix,1,multiply_by,Ptr_norm_factor) %>% t %>% Matrix
# # foo = Ptr_PtrEgr_norm_UMI_matrix %>%
# #     rbind(.,Ptr_PtrTar_norm_UMI_matrix[setdiff(rownames(Ptr_PtrTar_norm_UMI_matrix),
# #                                           rownames(.)),]) %>%
# #     rbind(.,Ptr_PtrLch_norm_UMI_matrix[setdiff(rownames(Ptr_PtrLch_norm_UMI_matrix),
# #                                           rownames(.)),])
# # stopifnot(all(foo == Ptr_PtrUnion3S_norm_UMI_matrix))
# 
# ##Input Egr data
# Egr_projection_umap = read.csv('Egr_projection.csv')
# selected_barcodes = Egr_projection_umap$Barcode
# 
# Egr_UMI_matrix = readRDS('RDS_Egr_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Egr_4S_UMI_matrix = readRDS('ortho_4S/RDS_Egr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Egr_PtrEgr_UMI_matrix = readRDS('ortho_Ptr_Egr/RDS_Egr_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# 
# Egr_norm_factor = 1000/apply(Egr_UMI_matrix,2,sum)
# Egr_norm_UMI_matrix = apply(Egr_UMI_matrix,1,multiply_by,Egr_norm_factor) %>% t %>% Matrix
# Egr_4S_norm_UMI_matrix = apply(Egr_4S_UMI_matrix,1,multiply_by,Egr_norm_factor) %>% t %>% Matrix
# Egr_PtrEgr_norm_UMI_matrix = apply(Egr_PtrEgr_UMI_matrix,1,multiply_by,Egr_norm_factor) %>% t %>% Matrix
# 
# ##Input Tar data
# Tar_projection_umap = read.csv('Tar_projection.csv')
# Tar_projection_umap$Barcode %<>% sub('Kun','Ta',.)
# selected_barcodes = Tar_projection_umap$Barcode
# 
# Tar_UMI_matrix = readRDS('RDS_Tar_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Tar_4S_UMI_matrix = readRDS('ortho_4S/RDS_Tar_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Tar_PtrTar_UMI_matrix = readRDS('ortho_Ptr_Tar/RDS_Tar_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# 
# Tar_norm_factor = 1000/apply(Tar_UMI_matrix,2,sum)
# Tar_norm_UMI_matrix = apply(Tar_UMI_matrix,1,multiply_by,Tar_norm_factor) %>% t %>% Matrix
# Tar_4S_norm_UMI_matrix = apply(Tar_4S_UMI_matrix,1,multiply_by,Tar_norm_factor) %>% t %>% Matrix
# Tar_PtrTar_norm_UMI_matrix = apply(Tar_PtrTar_UMI_matrix,1,multiply_by,Tar_norm_factor) %>% t %>% Matrix
# 
# ##Input Lch data
# Lch_projection_umap = readRDS('RDS_Lch_projection_umap.rds')
# selected_barcodes = Lch_projection_umap$Barcode
# 
# Lch_UMI_matrix = readRDS('RDS_Lch_SCseq_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Lch_4S_UMI_matrix = readRDS('ortho_4S/RDS_Lch_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# Lch_PtrLch_UMI_matrix = readRDS('ortho_Ptr_Lch/RDS_Lch_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes] %>% Matrix
# 
# Lch_norm_factor = 1000/apply(Lch_UMI_matrix,2,sum)
# Lch_norm_UMI_matrix = apply(Lch_UMI_matrix,1,multiply_by,Lch_norm_factor) %>% t %>% Matrix
# Lch_4S_norm_UMI_matrix = apply(Lch_4S_UMI_matrix,1,multiply_by,Lch_norm_factor) %>% t %>% Matrix
# Lch_PtrLch_norm_UMI_matrix = apply(Lch_PtrLch_UMI_matrix,1,multiply_by,Lch_norm_factor) %>% t %>% Matrix
# 
# save.image('Data_pack.RData')
load('Data_pack.RData')


#Define plotting function
plot_specific_feature = function(Species = 'Ptr',
                                 species_group = 'Ptr',
                                 Norm = T,
                                 point_size = 1,
                                 featureID = 'Potri.004G059600.v4.1',
                                 summary_title = 'Summary UMI',
                                 summary_method = 'none',
                                 output_figure = T,
                                 output_filename = paste0('Ptr_CellRanger_UMAP_',featureID,'.png'),
                                 color_scale = c(20,50),
                                 output_for_Jimmy = F){
    read_matrix_code = paste('plot_UMI_matrix =',Species)
    if(species_group %in% c('PtrEgr','PtrTar','PtrLch','4S','PtrUnion3S')){
        if(Species=='Ptr' & species_group%in%c('PtrEgr','PtrTar','PtrLch','PtrUnion3S')){
            read_matrix_code %<>% paste0('_PtrUnion3S')
        }else{
            read_matrix_code %<>% paste0('_',species_group)
        }
    }
    if(Norm)
        read_matrix_code %<>% paste0('_norm')
    read_matrix_code %<>% paste0('_UMI_matrix')
    eval(parse(text=read_matrix_code))
    
    read_projection_code = paste0('projection_umap =',Species,'_projection_umap')
    eval(parse(text=read_projection_code))
    
    pp1 = color_scale[1]
    pp2 = color_scale[2]
    COL = c(colorRampPalette(c('#C9D7EF','#eeeeec'))(pp1),
            colorRampPalette(c('#eeeeec','#ea4335'))(pp2-pp1),
            colorRampPalette(c('#ea4335','#ba1306'))(100-pp2))
    # Original blue: #4285f4
    
    if(output_figure) png(output_filename,
                          pointsize=10,width=20,height=15,units='cm',res=300)
    if(featureID %in% rownames(plot_UMI_matrix)){
        if(output_for_Jimmy) par(mai=c(0,0,0,0)) else par(mai=oriPar$mai)
        if(summary_method=='none'){
            plot_values = plot_UMI_matrix[featureID,] %>%
                as.numeric %>% add(1) %>% log2
        }else if(summary_method=='sum'){
            plot_values = plot_UMI_matrix[featureID,,drop=F] %>%
                apply(2,sum) %>% as.numeric %>% add(1) %>% log2
        }else if(summary_method=='max'){
            plot_values = plot_UMI_matrix[featureID,,drop=F] %>%
                apply(2,max) %>% as.numeric %>% add(1) %>% log2
        }
        COL_ind = round(plot_values/max(plot_values)*99+1)
        order_ind = order(plot_values)
        plot(type='n',
             projection_umap$UMAP.1,
             projection_umap$UMAP.2,
             main=ifelse(output_for_Jimmy,'',
                         paste0(species_group,': ',Species,'\n',
                                ifelse(length(featureID)>1,summary_title,featureID))),
             xlab='UMAP_1',ylab='UMAP_2',
             axes=!output_for_Jimmy)
        points(projection_umap$UMAP.1[order_ind],
               projection_umap$UMAP.2[order_ind],
               col=COL[COL_ind][order_ind],
               pch=20,cex=point_size)
        xlim = par('usr')[1:2]
        ylim = par('usr')[3:4]
        for(i in seq_along(COL)){
            x = seq(qunif(0.65,xlim[1],xlim[2]),
                    qunif(0.95,xlim[1],xlim[2]),
                    length.out=length(COL)+1)[c(i,i+1)]
            y = rep(ylim[2]+(ylim[2]-ylim[1])*7/100,2)
            lines(x,y,lwd=15,lend='butt',col=COL[i],xpd=T)
            if(i==pp1) lines(x,y,lwd=15,lend='butt',xpd=T)
            if(i==pp2) lines(x,y,lwd=15,lend='butt',xpd=T)
        }
        if(Norm){
            text(qunif(0.80,xlim[1],xlim[2]),
                 ylim[2]+(ylim[2]-ylim[1])*2.5/100,
                 'log2((UMI/totalUMI)*1000+1)',cex=1,xpd=T)
        }else{
            text(qunif(0.80,xlim[1],xlim[2]),
                 ylim[2]+(ylim[2]-ylim[1])*2.5/100,
                 'log2(UMI+1)',cex=1,xpd=T)
        }
        text(qunif(0.95,xlim[1],xlim[2]),
             ylim[2]+(ylim[2]-ylim[1])*12/100,
             sprintf('%.2f',max(plot_values)),cex=1.5,xpd=T)
        text(qunif(0.65,xlim[1],xlim[2]),
             ylim[2]+(ylim[2]-ylim[1])*12/100,
             sprintf('%.2f',0),cex=1.5,xpd=T)
        par(mai=oriPar$mai)
    }else{
        plot.new()
    }
    if(output_figure) dev.off()
}


plot_specific_feature_4in1 = function(Norm = T,
                                      point_size = 1,
                                      featureID = 'Cluster_0',
                                      summary_title = 'Summary UMI',
                                      summary_method = 'none',
                                      output_filename = paste0('4S_CellRanger_UMAP_',featureID,'.png'),
                                      color_scale = list(Ptr=c(20,50),Egr=c(20,50),Tar=c(20,50),Lch=c(20,50)),
                                      output_for_Jimmy = F){
    png(output_filename,
        pointsize=10,width=80,height=15,units='cm',res=300)
    par(mfcol=c(1,4))
    plot_specific_feature(Species = 'Ptr',
                          species_group = 'PtrUnion3S',
                          Norm = T,
                          point_size = point_size,
                          featureID = featureID,
                          summary_title = summary_title,
                          summary_method = summary_method,
                          output_figure = F,
                          output_filename = NULL,
                          color_scale = color_scale[['Ptr']],
                          output_for_Jimmy = output_for_Jimmy)
    for(sp in c('Egr','Tar','Lch')){
        plot_specific_feature(Species = sp,
                              species_group = paste0('Ptr',sp),
                              Norm = T,
                              point_size = point_size,
                              featureID = featureID,
                              summary_title = summary_title,
                              summary_method = summary_method,
                              output_figure = F,
                              output_filename = NULL,
                              color_scale = color_scale[[sp]],
                              output_for_Jimmy = output_for_Jimmy)
    }
    dev.off()
    par(oriPar)
}

save.image('Plot_specific_feature.RData')

