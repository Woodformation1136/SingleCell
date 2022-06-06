library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\Ptr_Osa')

#Setup the Seurat Object
Ptr.data = Read10X(data.dir = 'Ptr_filtered_feature_bc_matrix')
Ptr = CreateSeuratObject(counts = Ptr.data, project = 'Ptr', min.cells = 3, min.features = 200)
Ptr

Osa2021R1.data = Read10X(data.dir = 'Osa_2021R1_filtered_feature_bc_matrix')
Osa2021R2.data = Read10X(data.dir = 'Osa_2021R2_filtered_feature_bc_matrix')
Osa2021R1 = CreateSeuratObject(counts = Osa2021R1.data, project = 'Osa2021R1', min.cells = 3, min.features = 200)
Osa2021R2 = CreateSeuratObject(counts = Osa2021R2.data, project = 'Osa2021R2', min.cells = 3, min.features = 200)
Osa2021R1
Osa2021R2

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Osa2021R1 = NormalizeData(Osa2021R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Osa2021R2 = NormalizeData(Osa2021R2, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Osa2021R1 = FindVariableFeatures(Osa2021R1, selection.method = 'vst', nfeatures = 2000)
Osa2021R2 = FindVariableFeatures(Osa2021R2, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Osa2021R = FindIntegrationAnchors(object.list = list(Ptr,Osa2021R1,Osa2021R2),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Osa2021R = IntegrateData(anchorset = integration_anchors_Osa2021R,
                                         preserve.order = T,
                                         sample.tree = rbind(c(-3,-2),
                                                             c(-1,1)))


#Run the standard workflow for visualization and clustering
Combined_object_Osa2021R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Osa2021R %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)

# Combined_object_Osa2021R@reductions$umap@cell.embeddings[,2] %<>% multiply_by(-1)
rotation_angle = -45 * (pi/180)
rotation_matrix = matrix(c(cos(rotation_angle),-sin(rotation_angle),
                           sin(rotation_angle),cos(rotation_angle)),nrow=2,byrow=T)
Combined_object_Osa2021R@reductions$umap@cell.embeddings %<>%
    multiply_by_matrix(rotation_matrix) %>%
    set_colnames(c('UMAP_1','UMAP_2'))

saveRDS(Combined_object_Osa2021R,'RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Osa2021R.rds')
Combined_object_Osa2021R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Osa2021R.rds')



DimPlot(Combined_object_Osa2021R, reduction = 'umap', group.by = 'orig.ident', pt.size = 0.5)
DimPlot(Combined_object_Osa2021R, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.5)



Ptr_SCseq_cell_cluster_10 = readRDS('RDS_Ptr_SCseq_cell_cluster_10.rds')
COL = c('#1F77B4','#8C564B','#FF7F0F','#2AA02A','#F8E71C',
        '#9467BD','#D62728','#E377C2','#9B9B9B','#4B4B4B')
# Ptr_projection_umap_CellRanger = readRDS('RDS_Ptr_projection_umap.rds')
# plot(Ptr_projection_umap_CellRanger$UMAP.1,
#      Ptr_projection_umap_CellRanger$UMAP.2,
#      pch=20,col=COL[Ptr_SCseq_cell_cluster_10$Cluster])
Ptr_barcode2color = data.frame(Barcode = Ptr_SCseq_cell_cluster_10$Barcode,
                               Color = COL[Ptr_SCseq_cell_cluster_10$Cluster])


projection_UMAP = Combined_object_Osa2021R@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$Species = projection_UMAP %>% rownames %>% substr(1,3)
projection_UMAP$Seurat_clusters = Combined_object_Osa2021R@meta.data$seurat_clusters %>% as.numeric

Ptr_projection_UMAP = filter(projection_UMAP,Species=='Ptr')
Ptr_projection_UMAP$Barcode = Ptr_projection_UMAP %>% rownames %>% substring(5)
Osa_projection_UMAP = filter(projection_UMAP,Species=='Osa')

png('Ptr_Osa2021R_Ptr_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.5,
       col=with(Ptr_barcode2color,
                Color[match(Ptr_projection_UMAP$Barcode,Barcode)]))
dev.off()


png('Ptr_Osa2021R_Ptr_Seurat_UMAP_black.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(Ptr_projection_UMAP$UMAP_1,
     Ptr_projection_UMAP$UMAP_2,
     pch=20,col='black',cex=0.5,
     las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
dev.off()

png('Ptr_Osa2021R_Osa_Seurat_UMAP_gold.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(Osa_projection_UMAP$UMAP_1,
     Osa_projection_UMAP$UMAP_2,
     pch=20,col='#C59739',cex=0.3,
     las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Osa')
dev.off()

png('Ptr_Osa2021R_Both_Seurat_UMAP.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrOsa')
points(Osa_projection_UMAP$UMAP_1,
       Osa_projection_UMAP$UMAP_2,
       pch=20,col='#C59739',cex=0.3)
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,col='black',cex=0.5)
dev.off()

png('Ptr_Osa2021R_Both_Seurat_UMAP_colorful.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrOsa')
points(Osa_projection_UMAP$UMAP_1,
       Osa_projection_UMAP$UMAP_2,
       pch=20,col='#C59739',cex=0.3)
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.5,
       col=with(Ptr_barcode2color,
                Color[match(Ptr_projection_UMAP$Barcode,Barcode)]))
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
png('Ptr_Osa2021R_Both_Seurat_UMAP_colored_by_Seurat_cluster_Cre30.png',
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




##Check batch effect in Osa2021R
barcodes = rownames(Combined_object_Osa2021R@meta.data)
group_id = barcodes %>% strsplit('_') %>% sapply(function(x)ifelse(length(x)==2,
                                                                   x[1],paste(x[1:2],collapse='_')))
# Combined_object_Osa2021R@meta.data$orig.group = group_id
# DimPlot(Combined_object_Osa2021R, reduction = 'umap', group.by = 'orig.group', pt.size = 0.5)

projection_UMAP = Combined_object_Osa2021R@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$group_id = group_id

Ptr_projection_UMAP = filter(projection_UMAP,group_id=='Ptr')
Osa_2021R1_projection_UMAP = filter(projection_UMAP,group_id=='Osa_2021R1')
Osa_2021R2_projection_UMAP = filter(projection_UMAP,group_id=='Osa_2021R2')

png('Ptr_Osa2021R_Both_Seurat_UMAP_colored_by_batch.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
points(Osa_2021R2_projection_UMAP$UMAP_1,
       Osa_2021R2_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='#00bfc4')
points(Osa_2021R1_projection_UMAP$UMAP_1,
       Osa_2021R1_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col='#f8766d')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.5,
       col='gray')
legend(-10,10,bty='n',pch=20,cex=1.2,
       legend=c('Ptr','Osa_2021R1','Osa_2021R2'),
       col=c('gray','#f8766d','#00bfc4'))
dev.off()


#Highlight the xylem cells according to previous research
##Osa2021R
xylem_barcode = scan('Osa_2021R1_Xylem_barcodes.txt',what='')

projection_UMAP = Combined_object_Osa2021R@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$Species = projection_UMAP %>% rownames %>% substr(1,3)

Ptr_projection_UMAP = filter(projection_UMAP,Species=='Ptr')
Osa_projection_UMAP = filter(projection_UMAP,Species=='Osa')

png('Ptr_Osa2021R_Both_Seurat_UMAP_Ath_xylem.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrAth')
points(Osa_projection_UMAP$UMAP_1,
       Osa_projection_UMAP$UMAP_2,
       pch=20,col='#C59739',cex=0.3)
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.5,
       col='gray')
points(Osa_projection_UMAP[xylem_barcode,]$UMAP_1,
       Osa_projection_UMAP[xylem_barcode,]$UMAP_2,
       pch=20,col='#df35db',cex=0.5)
dev.off()

