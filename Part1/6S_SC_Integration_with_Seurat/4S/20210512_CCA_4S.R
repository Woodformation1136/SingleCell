library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\4S')

# 紫色: '#9467BD'
# 黃色: '#F8E71C'
# 棕色: '#8C564B'
# 藍色: '#1F77B4'
# 綠色: '#2AA02A'
# 橘色: '#FF7F0F'
# 紅色: '#D62728'
# 灰色: '#9B9B9B'
# 粉色: '#E377C2'
# 青色: '#50E3C2'
# 怪綠: '#63EE9B'
# 金色: '#C59739'
# 黃綠: '#C9DD31'
# 天藍: '#AEE9F4'
# 樵藍: '#46CBE5'
# 貴紫: '#A9517B'

#Setup the Seurat Object
Ptr.data = Read10X(data.dir = 'Ptr_filtered_feature_bc_matrix')
Lch.data = Read10X(data.dir = 'Lch_filtered_feature_bc_matrix')
Egr.data = Read10X(data.dir = 'Egr_filtered_feature_bc_matrix')
Tar.data = Read10X(data.dir = 'Tar_filtered_feature_bc_matrix')
Ptr = CreateSeuratObject(counts = Ptr.data, project = 'Ptr', min.cells = 3, min.features = 200)
Lch = CreateSeuratObject(counts = Lch.data, project = 'Lch', min.cells = 3, min.features = 200)
Egr = CreateSeuratObject(counts = Egr.data, project = 'Egr', min.cells = 3, min.features = 200)
Tar = CreateSeuratObject(counts = Tar.data, project = 'Tar', min.cells = 3, min.features = 200)
# 9601 gene ortholog clusters in total
Ptr
Lch
Egr
Tar

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Lch = NormalizeData(Lch, normalization.method = 'LogNormalize', scale.factor = 10000)
Egr = NormalizeData(Egr, normalization.method = 'LogNormalize', scale.factor = 10000)
Tar = NormalizeData(Tar, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Lch = FindVariableFeatures(Lch, selection.method = 'vst', nfeatures = 2000)
Egr = FindVariableFeatures(Egr, selection.method = 'vst', nfeatures = 2000)
Tar = FindVariableFeatures(Tar, selection.method = 'vst', nfeatures = 2000)

#Select variable genes common to both samples
# Ptr_variables = VariableFeatures(Ptr)
# Lch_variables = VariableFeatures(Lch)
# genes.use = intersect(Ptr_variables,Lch_variables)

#Scale the data
# Ptr = ScaleData(Ptr, features = rownames(Ptr))
# Lch = ScaleData(Lch, features = rownames(Lch))

#Perform Canonical Correlation Analysis (CCA)
# Ptr_subset = subset(Ptr, cells = colnames(Ptr)[1:40])
# Lch_subset = subset(Lch, cells = colnames(Lch)[1:40])
# Combined_object = RunCCA(object1 = Ptr, object2 = Lch)
# str(Combined_object)
# DimPlot(object=Combined_object,
#         reduction='cca',
#         group.by='orig.ident')

#Perform l2 normalization on CCs
# Combined_object = L2CCA(Combined_object)
# DimPlot(object=Combined_object,
#         reduction='cca.l2',
#         group.by='orig.ident')

#Find integration anchors and integrate data
integration_anchors = FindIntegrationAnchors(object.list = list(Ptr,Egr,Tar,Lch),
                                             anchor.features = 2000,
                                             scale = TRUE,
                                             reduction = 'cca',
                                             l2.norm = TRUE,
                                             k.anchor = 5)

Combined_object = IntegrateData(anchorset = integration_anchors)

# saveRDS(integration_anchors,'RDS_integration_anchors_SF10000_AF2000_KA5.rds')
integration_anchors = readRDS('RDS_integration_anchors_SF10000_AF2000_KA5.rds')


#Run the standard workflow for visualization and clustering
# stopifnot(DefaultAssay(Combined_object) == 'integrated')
Combined_object = ScaleData(Combined_object)
Combined_object = RunPCA(Combined_object, npcs = 30)
Combined_object = RunUMAP(Combined_object, reduction = 'pca', dims = 1:30)
Combined_object = FindNeighbors(Combined_object, reduction = 'pca', dims = 1:30, k.param = 3)
Combined_object = FindClusters(Combined_object, resolution = 0.5)

# saveRDS(Combined_object,'RDS_Combined_object_SF10000_AF2000_KA5.rds')
Combined_object = readRDS('RDS_Combined_object_SF10000_AF2000_KA5.rds')


#Plot and save the Seurat UMAP
##Color with species
p1 = DimPlot(Combined_object, reduction = 'umap', group.by = 'orig.ident', pt.size = 0.5)
tiff('Four_species.tiff',
     width=20,height=15,units='cm',
     res=600)
p1
dev.off()

##color with Seurat cluster
p2 = DimPlot(Combined_object, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.5)
tiff('Four_species_clusters.tiff',
     width=40,height=15,units='cm',
     res=600)
p1 + p2
dev.off()

##color with specific feature
# DefaultAssay(Combined_object) = 'RNA'
# feature='Cluster-2985'
# tiff(paste0('Four_species_',feature,'.tiff'),
#      width=48,height=12,units='cm',
#      res=300)
# FeaturePlot(Combined_object,
#             features=feature,
#             pt.size=0.5,
#             order=T,
#             slot='data',
#             cols=brewer.pal(9,'YlOrRd'),
#             reduction='umap',
#             split.by='orig.ident')
# dev.off()


#Combime with CellRanger data
##Input the UMAP and clustering data
Ptr_SCseq_cell_cluster_10 = readRDS('RDS_Ptr_SCseq_cell_cluster_10.rds')
Lch_SCseq_cell_cluster_10 = readRDS('RDS_Lch_SCseq_cell_cluster_10.rds')
Egr_SCseq_cell_cluster_10 = read.csv('Egr_clusters.csv')
Tar_SCseq_cell_cluster_18 = read.csv('Tar_clusters_18.csv')
Tar_SCseq_cell_cluster_18$Barcode %<>% sub('Kun','Ta',.)
Ptr_projection_umap_CellRanger = readRDS('RDS_Ptr_projection_umap.rds')
Lch_projection_umap_CellRanger = readRDS('RDS_Lch_projection_umap.rds')
Egr_projection_umap_CellRanger = read.csv('Egr_projection.csv')
Tar_projection_umap_CellRanger = read.csv('Tar_projection.csv')
Tar_projection_umap_CellRanger$Barcode %<>% sub('Kun','Ta',.)

##Replot the CellRanger UMAP and color with CellRanger clusters and 
##build the barcode2color table
COL = c('#1F77B4','#8C564B','#FF7F0F','#2AA02A','#F8E71C',
        '#9467BD','#D62728','#E377C2','#9B9B9B','#4B4B4B')
# png('Ptr_CellRanger_UMAP.png',
#     pointsize=10,width=20,height=15,units='cm',res=300)
# par(mar=c(0,0,0,0))
plot(Ptr_projection_umap_CellRanger$UMAP.1,
     Ptr_projection_umap_CellRanger$UMAP.2,
     pch=20,cex=1.2,
     col=COL[Ptr_SCseq_cell_cluster_10$Cluster])
# par(oriPar)
# dev.off()
Ptr_barcode2color = data.frame(Barcode = Ptr_SCseq_cell_cluster_10$Barcode,
                               Color = COL[Ptr_SCseq_cell_cluster_10$Cluster])

COL = c('#9467BD','#F8E71C','#8C564B','#1F77B4','#2AA02A',
        '#FF7F0F','#D62728','#E377C2','#4B4B4B','#9B9B9B')
plot(Egr_projection_umap_CellRanger$UMAP.1,
     Egr_projection_umap_CellRanger$UMAP.2,
     pch=20,col=COL[Egr_SCseq_cell_cluster_10$Cluster])
Egr_barcode2color = data.frame(Barcode = Egr_SCseq_cell_cluster_10$Barcode,
                               Color = COL[Egr_SCseq_cell_cluster_10$Cluster])

COL = c('#1F77B4','#8C564B','#D62728','#F8E71C','#FF7F0F',
        '#9467BD','#E377C2','#2AA02A','#63EE9B','#50E3C2')
plot(Lch_projection_umap_CellRanger$UMAP.1,
     Lch_projection_umap_CellRanger$UMAP.2,
     pch=20,col=COL[Lch_SCseq_cell_cluster_10$Cluster])
Lch_barcode2color = data.frame(Barcode = Lch_SCseq_cell_cluster_10$Barcode,
                               Color = COL[Lch_SCseq_cell_cluster_10$Cluster])

COL = c('#1F77B4','#E377C2','#55A3FF','#8C564B','#9467BD',
        '#FF7F0F','#46CBE5','#F8E71C','#2AA02A','#9B9B9B',
        '#9B9B9B','#9B9B9B','#9B9B9B','#9B9B9B','#9B9B9B',
        '#9B9B9B','#9B9B9B','#9B9B9B')
plot(Tar_projection_umap_CellRanger$UMAP.1,
     Tar_projection_umap_CellRanger$UMAP.2,
     pch=20,col=COL[Tar_SCseq_cell_cluster_18$Cluster])
Tar_barcode2color = data.frame(Barcode = Tar_SCseq_cell_cluster_18$Barcode,
                               Color = COL[Tar_SCseq_cell_cluster_18$Cluster])



##Prepare the table to plot the Seurat UMAP
projection_UMAP = Combined_object@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$SpBarcode = projection_UMAP %>% rownames
projection_UMAP$Species = projection_UMAP %>% rownames %>% substr(1,3)
projection_UMAP$Barcode = projection_UMAP %>% rownames %>% substring(5)
projection_UMAP$Seurat_clusters = Combined_object@meta.data$seurat_clusters %>% as.numeric

total_barcode2color = data.frame(SpBarcode = c(paste0('Ptr_',Ptr_barcode2color$Barcode),
                                               paste0('Egr_',Egr_barcode2color$Barcode),
                                               paste0('Tar_',Tar_barcode2color$Barcode),
                                               paste0('Lch_',Lch_barcode2color$Barcode)),
                                 CellRanger_color = c(Ptr_barcode2color$Color,
                                                      Egr_barcode2color$Color,
                                                      Tar_barcode2color$Color,
                                                      Lch_barcode2color$Color))
projection_UMAP %<>% left_join(total_barcode2color, by = 'SpBarcode') %>%
    set_rownames(projection_UMAP$SpBarcode)
stopifnot(all(!is.na(projection_UMAP$CellRanger_color)))

Ptr_projection_UMAP = filter(projection_UMAP,Species=='Ptr')
Egr_projection_UMAP = filter(projection_UMAP,Species=='Egr')
Tar_projection_UMAP = filter(projection_UMAP,Species=='Tar')
Lch_projection_UMAP = filter(projection_UMAP,Species=='Lch')

##Plot the Seurat UMAP color with CellRanger cluster
png('4S_Ptr_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=8,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=1,
       col=Ptr_projection_UMAP$CellRanger_color)
dev.off()

png('4S_Egr_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=8,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Egr')
points(Egr_projection_UMAP$UMAP_1,
       Egr_projection_UMAP$UMAP_2,
       pch=20,cex=1,
       col=Egr_projection_UMAP$CellRanger_color)
dev.off()

png('4S_Tar_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=8,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Tar')
points(Tar_projection_UMAP$UMAP_1,
       Tar_projection_UMAP$UMAP_2,
       pch=20,cex=1,
       col=Tar_projection_UMAP$CellRanger_color)
dev.off()

png('4S_Lch_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=8,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Lch')
points(Lch_projection_UMAP$UMAP_1,
       Lch_projection_UMAP$UMAP_2,
       pch=20,cex=1,
       col=Lch_projection_UMAP$CellRanger_color)
dev.off()


#Try to plot the Seurat UMAP with Seurat clusters matching the CellRanger colors
##Construct the color table based on the major CellRanger color within each cluster
foo = with(projection_UMAP,table(Seurat_clusters,CellRanger_color))
foo_dimnames = attr(foo,'dimnames')
color_table = data.frame(Seurat_cluster = as.numeric(foo_dimnames$Seurat_clusters),
                         Color = foo_dimnames$CellRanger_color[apply(foo,1,which.max)])
Seurat_cluster_center = sapply(color_table$Seurat_cluster,
                               function(i){
                                   foo = filter(projection_UMAP,Seurat_clusters==i)
                                   return(c(mean(foo$UMAP_1),mean(foo$UMAP_2)))
                               })

png('4S_trying_colors.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(projection_UMAP$UMAP_1,
     projection_UMAP$UMAP_2,
     pch=20,cex=0.5,
     col=with(color_table,
              Color[match(projection_UMAP$Seurat_clusters,Seurat_cluster)]),
     xlab='UMAP_1',ylab='UMAP_2',main='4S')
legend(9,6.5,xpd=T,bty='n',
       pch=20,col=color_table$Color,
       legend=color_table$Seurat_cluster,pt.cex=2)
for(i in 1:nrow(color_table)) text(Seurat_cluster_center[1,i],
                                   Seurat_cluster_center[2,i],
                                   color_table$Seurat_cluster[i],cex=1.5)
dev.off()


############# TO BE CONTINUED ####################


##Change some color by hand
handcraft_color_table = color_table
handcraft_color_table[c(14,17,20,23,24,26,29),'Color'] = '#9B9B9B'

png('PtrEgr_Both_Seurat_UMAP_colored_by_Seurat_cluster.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(projection_UMAP$UMAP_1,
     projection_UMAP$UMAP_2,
     pch=20,cex=0.3,
     col=with(handcraft_color_table,
              Color[match(projection_UMAP$Seurat_clusters,Seurat_cluster)]),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrEgr')
for(i in 1:nrow(color_table)) text(Seurat_cluster_center[1,i],
                                   Seurat_cluster_center[2,i],
                                   color_table$Seurat_cluster[i],cex=1.5)
dev.off()

##Color with the Cre-picking color
Cre30_color_table = data.frame(Seurat_cluster = 1:30,
                               Color = c('#BDC2C5','#DF35DB','#FFFF55','#FE3BAC','#D6DCE5',
                                         '#4F3098','#FF71AD','#4F9BF9','#148F7C','#951DA2',
                                         '#A74831','#A262F0','#5C4C32','#ADEB4D','#FFE21C',
                                         '#FD3A5C','#2695EA','#87ABFD','#FE4496','#44BF77',
                                         '#F0AB2C','#FF3E3F','#FE7D1C','#16EE91','#6DC6FD',
                                         '#69D8E9','#DF2036','#BB986A','#282DE8','#595959'))
png('PtrEgr_Both_Seurat_UMAP_colored_by_Seurat_cluster_Cre30.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(projection_UMAP$UMAP_1,
     projection_UMAP$UMAP_2,
     pch=20,cex=0.3,
     col=with(Cre30_color_table,
              Color[match(projection_UMAP$Seurat_clusters,Seurat_cluster)]),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrEgr')
for(i in 1:nrow(color_table)) text(Seurat_cluster_center[1,i],
                                   Seurat_cluster_center[2,i],
                                   color_table$Seurat_cluster[i],cex=1.5)
dev.off()


#Conduct the differential expression analysis (DEA) for each handcraft-color-cluster
projection_UMAP$Seurat_color = with(handcraft_color_table,
                                    Color[match(as.numeric(Combined_object@meta.data$seurat_clusters),
                                                Seurat_cluster)])
saveRDS(projection_UMAP,'RDS_4S_Seurat_projection_UMAP.rds')
dim(projection_UMAP) # 9556 8

Combined_object = SetIdent(Combined_object,
                           value=projection_UMAP$Seurat_color)
##One-versus-all DEA
all.markers = FindAllMarkers(object=Combined_object, test.use='MAST', only.pos=T)
all.markers %>% str # 2598 7
saveRDS(all.markers,'RDS_PtrEgr_allmarkers.rds')

##Take two for checks
filter(all.markers,cluster=='#D62728') %>% head
FeaturePlot(Combined_object,features='Cluster-2455',pt.size=1)
FeaturePlot(Combined_object,features='Cluster-7685',pt.size=1)

