library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\Ptr_Ath')

#Setup the Seurat Object
Ptr.data = Read10X(data.dir = 'Ptr_filtered_feature_bc_matrix')
Ptr = CreateSeuratObject(counts = Ptr.data, project = 'Ptr', min.cells = 3, min.features = 200)
Ptr

Ath2021L.data = Read10X(data.dir = 'Ath_2021L_filtered_feature_bc_matrix')
Ath2021L = CreateSeuratObject(counts = Ath2021L.data, project = 'Ath2021L', min.cells = 3, min.features = 200)
Ath2021L

Ath2019R1.data = Read10X(data.dir = 'Ath_2019R1_filtered_feature_bc_matrix')
Ath2019R2.data = Read10X(data.dir = 'Ath_2019R2_filtered_feature_bc_matrix')
Ath2019R3.data = Read10X(data.dir = 'Ath_2019R3_filtered_feature_bc_matrix')
Ath2019R1 = CreateSeuratObject(counts = Ath2019R1.data, project = 'Ath2019R1', min.cells = 3, min.features = 200)
Ath2019R2 = CreateSeuratObject(counts = Ath2019R2.data, project = 'Ath2019R2', min.cells = 3, min.features = 200)
Ath2019R3 = CreateSeuratObject(counts = Ath2019R3.data, project = 'Ath2019R3', min.cells = 3, min.features = 200)
Ath2019R1
Ath2019R2
Ath2019R3

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2021L = NormalizeData(Ath2021L, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R1 = NormalizeData(Ath2019R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R2 = NormalizeData(Ath2019R2, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R3 = NormalizeData(Ath2019R3, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Ath2021L = FindVariableFeatures(Ath2021L, selection.method = 'vst', nfeatures = 2000)
Ath2019R1 = FindVariableFeatures(Ath2019R1, selection.method = 'vst', nfeatures = 2000)
Ath2019R2 = FindVariableFeatures(Ath2019R2, selection.method = 'vst', nfeatures = 2000)
Ath2019R3 = FindVariableFeatures(Ath2019R3, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Ath2021L = FindIntegrationAnchors(object.list = list(Ptr,Ath2021L),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)
integration_anchors_Ath2019R = FindIntegrationAnchors(object.list = list(Ptr,Ath2019R1,Ath2019R2,Ath2019R3),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Ath2021L = IntegrateData(anchorset = integration_anchors_Ath2021L)
Combined_object_Ath2019R = IntegrateData(anchorset = integration_anchors_Ath2019R,
                                         preserve.order = T,
                                         sample.tree = rbind(c(-1,-2),
                                                             c(1,-3),
                                                             c(2,-4)))


#Run the standard workflow for visualization and clustering
Combined_object_Ath2021L %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Ath2021L %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)
Combined_object_Ath2021L@reductions$umap@cell.embeddings[,2] %<>% multiply_by(-1)
rotation_angle = 60 * (pi/180)
rotation_matrix = matrix(c(cos(rotation_angle),-sin(rotation_angle),
                           sin(rotation_angle),cos(rotation_angle)),nrow=2,byrow=T)
Combined_object_Ath2021L@reductions$umap@cell.embeddings %<>% 
    multiply_by_matrix(rotation_matrix) %>% 
    set_colnames(c('UMAP_1','UMAP_2'))

saveRDS(Combined_object_Ath2021L,'RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2021L.rds')
Combined_object_Ath2021L = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2021L.rds')


Combined_object_Ath2019R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Ath2019R %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)
Combined_object_Ath2019R@reductions$umap@cell.embeddings[,1] %<>% multiply_by(-1)

saveRDS(Combined_object_Ath2019R,'RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2019R.rds')
Combined_object_Ath2019R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2019R.rds')


Ath_output_various_integrated_umap = function(Combined_object = Combined_object_Ath2021L,
                                              plot_prefix = 'Ptr_Ath2021L_'){
    # p1 = DimPlot(Combined_object, reduction = 'umap', group.by = 'orig.ident', pt.size = 0.5)
    # p2 = DimPlot(Combined_object, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.5)
    # 
    # tiff(paste0(plot_prefix,'Both.tiff'),
    #      width=20,height=15,units='cm',
    #      res=600)
    # p1
    # dev.off()
    # 
    # tiff(paste0(plot_prefix,'clusters.tiff'),
    #      width=40,height=15,units='cm',
    #      res=600)
    # p1 + p2
    # dev.off()
    # 
    # tiff(paste0(plot_prefix,'Both_Seurat_UMAP_colored_by_Seurat_cluster.tiff'),
    #      width=20,height=15,units='cm',
    #      res=300,pointsize=10)
    # p2
    # dev.off()
    
    
    Ptr_SCseq_cell_cluster_10 = readRDS('RDS_Ptr_SCseq_cell_cluster_10.rds')
    Ptr_projection_umap_CellRanger = readRDS('RDS_Ptr_projection_umap.rds')
    COL = c('#1F77B4','#8C564B','#FF7F0F','#2AA02A','#F8E71C',
            '#9467BD','#D62728','#E377C2','#9B9B9B','#4B4B4B')
    # plot(Ptr_projection_umap_CellRanger$UMAP.1,
    #      Ptr_projection_umap_CellRanger$UMAP.2,
    #      pch=20,col=COL[Ptr_SCseq_cell_cluster_10$Cluster])
    Ptr_barcode2color = data.frame(Barcode = Ptr_SCseq_cell_cluster_10$Barcode,
                                   Color = COL[Ptr_SCseq_cell_cluster_10$Cluster])
    
    
    projection_UMAP = Combined_object@reductions$umap@cell.embeddings %>% as.data.frame
    projection_UMAP$Species = projection_UMAP %>% rownames %>% substr(1,3)
    projection_UMAP$Barcode = projection_UMAP %>% rownames %>% substring(5)
    projection_UMAP$Seurat_clusters = Combined_object@meta.data$seurat_clusters %>% as.numeric
    
    Ptr_projection_UMAP = filter(projection_UMAP,Species=='Ptr')
    Ath_projection_UMAP = filter(projection_UMAP,Species=='Ath')
    
    png(paste0(plot_prefix,'Ptr_Seurat_UMAP_colored_by_CellRanger_cluster.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot(0,0,type='n',las=1,
         xlim=range(projection_UMAP$UMAP_1),
         ylim=range(projection_UMAP$UMAP_2),
         xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
    # abline(h=seq(-5,15,5),col='gray')
    # abline(v=seq(-10,5,5),col='gray')
    points(Ptr_projection_UMAP$UMAP_1,
           Ptr_projection_UMAP$UMAP_2,
           pch=20,cex=0.3,
           col=with(Ptr_barcode2color,
                    Color[match(Ptr_projection_UMAP$Barcode,Barcode)]))
    dev.off()
    
    
    png(paste0(plot_prefix,'Ptr_Seurat_UMAP_black.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot(Ptr_projection_UMAP$UMAP_1,
         Ptr_projection_UMAP$UMAP_2,
         pch=20,col='black',cex=0.3,
         las=1,
         xlim=range(projection_UMAP$UMAP_1),
         ylim=range(projection_UMAP$UMAP_2),
         xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
    dev.off()
    
    png(paste0(plot_prefix,'Ath_Seurat_UMAP_gold.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot(Ath_projection_UMAP$UMAP_1,
         Ath_projection_UMAP$UMAP_2,
         pch=20,col='#C59739',cex=0.3,
         las=1,
         xlim=range(projection_UMAP$UMAP_1),
         ylim=range(projection_UMAP$UMAP_2),
         xlab='UMAP_1',ylab='UMAP_2',main='Ath')
    dev.off()
    
    png(paste0(plot_prefix,'Both_Seurat_UMAP.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot(0,0,type='n',las=1,
         xlim=range(projection_UMAP$UMAP_1),
         ylim=range(projection_UMAP$UMAP_2),
         xlab='UMAP_1',ylab='UMAP_2',main='PtrAth')
    points(Ptr_projection_UMAP$UMAP_1,
           Ptr_projection_UMAP$UMAP_2,
           pch=20,col='black',cex=0.3)
    points(Ath_projection_UMAP$UMAP_1,
           Ath_projection_UMAP$UMAP_2,
           pch=20,col='#C59739',cex=0.3)
    dev.off()
    
    png(paste0(plot_prefix,'Both_Seurat_UMAP_colorful.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot(0,0,type='n',las=1,
         xlim=range(projection_UMAP$UMAP_1),
         ylim=range(projection_UMAP$UMAP_2),
         xlab='UMAP_1',ylab='UMAP_2',main='PtrAth')
    points(Ptr_projection_UMAP$UMAP_1,
           Ptr_projection_UMAP$UMAP_2,
           pch=20,cex=0.3,
           col=with(Ptr_barcode2color,
                    Color[match(Ptr_projection_UMAP$Barcode,Barcode)]))
    points(Ath_projection_UMAP$UMAP_1,
           Ath_projection_UMAP$UMAP_2,
           pch=20,col='#C59739',cex=0.3)
    dev.off()
    
    ##Color with the Cre-picking color
    n_cluster = max(projection_UMAP$Seurat_clusters)
    Seurat_cluster_center = sapply(seq(n_cluster),
                                   function(i){
                                       foo = filter(projection_UMAP,Seurat_clusters==i)
                                       return(c(mean(foo$UMAP_1),mean(foo$UMAP_2)))
                                   })
    Cre30_color_table = data.frame(Seurat_cluster = seq(n_cluster),
                                   Color = scales::hue_pal()(n_cluster))
    png(paste0(plot_prefix,'Both_Seurat_UMAP_colored_by_Seurat_cluster_Cre30.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot(projection_UMAP$UMAP_1,
         projection_UMAP$UMAP_2,
         pch=20,cex=0.3,
         col=with(Cre30_color_table,
                  Color[match(projection_UMAP$Seurat_clusters,Seurat_cluster)]),
         xlab='UMAP_1',ylab='UMAP_2',main='PtrAth')
    for(i in seq(n_cluster)) text(Seurat_cluster_center[1,i],
                                  Seurat_cluster_center[2,i],
                                  i,cex=2.5)
    dev.off()
}


Ath_output_various_integrated_umap(Combined_object = Combined_object_Ath2021L,
                                   plot_prefix = 'Ptr_Ath2021L_')
Ath_output_various_integrated_umap(Combined_object = Combined_object_Ath2019R,
                                   plot_prefix = 'Ptr_Ath2019R_')


##Check batch effect in Ath_2019R
barcodes = rownames(Combined_object_Ath2019R@meta.data)
group_id = barcodes %>% strsplit('_') %>% sapply(function(x)ifelse(length(x)==2,
                                                                   x[1],paste(x[1:2],collapse='_')))
Combined_object_Ath2019R@meta.data$orig.group = group_id
DimPlot(Combined_object_Ath2019R, reduction = 'umap', group.by = 'orig.group', pt.size = 0.5)

projection_UMAP = Combined_object_Ath2019R@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$group_id = group_id

Ptr_projection_UMAP = filter(projection_UMAP,group_id=='Ptr')
Ath_2019R1_projection_UMAP = filter(projection_UMAP,group_id=='Ath_2019R1')
Ath_2019R2_projection_UMAP = filter(projection_UMAP,group_id=='Ath_2019R2')
Ath_2019R3_projection_UMAP = filter(projection_UMAP,group_id=='Ath_2019R3')

png('Ptr_Ath2019R_Both_Seurat_UMAP_colored_by_batch.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='gray')
points(Ath_2019R1_projection_UMAP$UMAP_1,
       Ath_2019R1_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='#f8766d')
points(Ath_2019R2_projection_UMAP$UMAP_1,
       Ath_2019R2_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='#7cae00')
points(Ath_2019R3_projection_UMAP$UMAP_1,
       Ath_2019R3_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='#00bfc4')
legend(6,-8,bty='n',pch=20,cex=1.2,
       legend=c('Ptr','Ath_2019R1','Ath_2019R2','Ath_2019R3'),
       col=c('gray','#f8766d','#7cae00','#00bfc4'))
dev.off()


#Highlight the xylem cells according to previous research
##Ath2021L
cluster_table = read.csv('Ptr_Ath2021L_K-Means 20.csv')
xylem_barcode = with(cluster_table,Barcode[X20=='Cluster 17'])

projection_UMAP = Combined_object_Ath2021L@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$Species = projection_UMAP %>% rownames %>% substr(1,3)

Ptr_projection_UMAP = filter(projection_UMAP,Species=='Ptr')
Ath_projection_UMAP = filter(projection_UMAP,Species=='Ath')

png('Ptr_Ath2021L_Both_Seurat_UMAP_Ath_xylem.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrAth')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='gray')
points(Ath_projection_UMAP$UMAP_1,
       Ath_projection_UMAP$UMAP_2,
       pch=20,col='#C59739',cex=0.3)
points(Ath_projection_UMAP[paste0('Ath_2021L_',xylem_barcode),]$UMAP_1,
       Ath_projection_UMAP[paste0('Ath_2021L_',xylem_barcode),]$UMAP_2,
       pch=20,col='#df35db',cex=0.6)
dev.off()

##Ath2019R
xylem_barcode = scan('Ptr_Ath2019R_Xylem_barcodes.txt',what='')

projection_UMAP = Combined_object_Ath2019R@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$Species = projection_UMAP %>% rownames %>% substr(1,3)

Ptr_projection_UMAP = filter(projection_UMAP,Species=='Ptr')
Ath_projection_UMAP = filter(projection_UMAP,Species=='Ath')

png('Ptr_Ath2019R_Both_Seurat_UMAP_Ath_xylem.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrAth')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='gray')
points(Ath_projection_UMAP$UMAP_1,
       Ath_projection_UMAP$UMAP_2,
       pch=20,col='#C59739',cex=0.3)
points(Ath_projection_UMAP[xylem_barcode,]$UMAP_1,
       Ath_projection_UMAP[xylem_barcode,]$UMAP_2,
       pch=20,col='#df35db',cex=0.6)
dev.off()

