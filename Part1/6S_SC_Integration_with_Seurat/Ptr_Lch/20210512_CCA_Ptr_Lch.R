library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\Single_cell\\20210600_SCseq_species_integration\\20210600_Ptr_Lch')

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
saveRDS(Ptr.data,'RDS_PtrLch_Ptr_SCseq_ortho_UMI_matrix.rds')
saveRDS(Lch.data,'RDS_PtrLch_Lch_SCseq_ortho_UMI_matrix.rds')

Ptr = CreateSeuratObject(counts = Ptr.data, project = 'Ptr', min.cells = 3, min.features = 200)
Lch = CreateSeuratObject(counts = Lch.data, project = 'Lch', min.cells = 3, min.features = 200)
Ptr
Lch

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Lch = NormalizeData(Lch, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Lch = FindVariableFeatures(Lch, selection.method = 'vst', nfeatures = 2000)

#Find integration anchors and integrate data
integration_anchors = FindIntegrationAnchors(object.list = list(Ptr,Lch),
                                             anchor.features = 2000,
                                             scale = TRUE,
                                             reduction = 'cca',
                                             l2.norm = TRUE,
                                             k.anchor = 5)

Combined_object = IntegrateData(anchorset = integration_anchors)

# saveRDS(integration_anchors,'RDS_integration_anchors_SF10000_AF2000_KA5_Ptr_Lch.rds')
integration_anchors = readRDS('RDS_integration_anchors_SF10000_AF2000_KA5_Ptr_Lch.rds')


#Run the standard workflow for visualization and clustering
# stopifnot(DefaultAssay(Combined_object) == 'integrated')
Combined_object = ScaleData(Combined_object)
Combined_object = RunPCA(Combined_object, npcs = 30)
Combined_object = RunUMAP(Combined_object, reduction = 'pca', dims = 1:30)
Combined_object = FindNeighbors(Combined_object, reduction = 'pca', dims = 1:30, k.param = 3)
Combined_object = FindClusters(Combined_object, resolution = 0.5)

# saveRDS(Combined_object,'RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Lch.rds')
Combined_object = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Lch.rds')


#Plot and save the Seurat UMAP
##Color with species
p1 = DimPlot(Combined_object, reduction = 'umap', group.by = 'orig.ident', pt.size = 0.5)
tiff('Ptr_Lch.tiff',
     width=20,height=15,units='cm',
     res=600)
p1
dev.off()

##color with Seurat cluster
p2 = DimPlot(Combined_object, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.5)
tiff('Ptr_Lch_clusters.tiff',
     width=40,height=15,units='cm',
     res=600)
p1 + p2
dev.off()

##color with specific feature
# DefaultAssay(Combined_object) = 'RNA'
# feature = 'Cluster-5814'
# tiff(paste0('PtrLch_',feature,'.tiff'),
#      width=32,height=12,units='cm',
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
Ptr_projection_umap_CellRanger = readRDS('RDS_Ptr_projection_umap.rds')
Lch_projection_umap_CellRanger = readRDS('RDS_Lch_projection_umap.rds')

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

COL = c('#1F77B4','#8C564B','#D62728','#F8E71C','#FF7F0F',
        '#9467BD','#E377C2','#2AA02A','#63EE9B','#50E3C2')
plot(Lch_projection_umap_CellRanger$UMAP.1,
     Lch_projection_umap_CellRanger$UMAP.2,
     pch=20,col=COL[Lch_SCseq_cell_cluster_10$Cluster])
Lch_barcode2color = data.frame(Barcode = Lch_SCseq_cell_cluster_10$Barcode,
                               Color = COL[Lch_SCseq_cell_cluster_10$Cluster])

##Prepare the table to plot the Seurat UMAP
projection_UMAP = Combined_object@reductions$umap@cell.embeddings %>% as.data.frame
projection_UMAP$SpBarcode = projection_UMAP %>% rownames
projection_UMAP$Species = projection_UMAP %>% rownames %>% substr(1,3)
projection_UMAP$Barcode = projection_UMAP %>% rownames %>% substring(5)
projection_UMAP$Seurat_clusters = Combined_object@meta.data$seurat_clusters %>% as.numeric

total_barcode2color = data.frame(SpBarcode = c(paste0('Ptr_',Ptr_barcode2color$Barcode),
                                               paste0('Lch_',Lch_barcode2color$Barcode)),
                                 CellRanger_color = c(Ptr_barcode2color$Color,
                                                      Lch_barcode2color$Color))
projection_UMAP %<>% left_join(total_barcode2color, by = 'SpBarcode') %>%
    set_rownames(projection_UMAP$SpBarcode)
stopifnot(all(!is.na(projection_UMAP$CellRanger_color)))

Ptr_projection_UMAP = filter(projection_UMAP,Species=='Ptr')
Lch_projection_UMAP = filter(projection_UMAP,Species=='Lch')

##Plot the Seurat UMAP color with CellRanger cluster
png('PtrLch_Ptr_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col=Ptr_projection_UMAP$CellRanger_color)
dev.off()

png('PtrLch_Lch_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Lch')
points(Lch_projection_UMAP$UMAP_1,
       Lch_projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col=Lch_projection_UMAP$CellRanger_color)
dev.off()

png('PtrLch_Both_Seurat_UMAP_colored_by_CellRanger_cluster.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrLch')
points(projection_UMAP$UMAP_1,
       projection_UMAP$UMAP_2,
       pch=20,cex=0.3,
       col=projection_UMAP$CellRanger_color)
dev.off()

##Plot the Seurat UMAP color with black and gold
png('PtrLch_Ptr_Seurat_UMAP_black.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(Ptr_projection_UMAP$UMAP_1,
     Ptr_projection_UMAP$UMAP_2,
     pch=20,col='black',cex=0.3,
     las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Ptr')
dev.off()

png('PtrLch_Lch_Seurat_UMAP_gold.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(Lch_projection_UMAP$UMAP_1,
     Lch_projection_UMAP$UMAP_2,
     pch=20,col='#C59739',cex=0.3,
     las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='Lch')
dev.off()

png('PtrLch_Both_Seurat_UMAP.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrLch')
points(Ptr_projection_UMAP$UMAP_1,
       Ptr_projection_UMAP$UMAP_2,
       pch=20,col='black',cex=0.3)
points(Lch_projection_UMAP$UMAP_1,
       Lch_projection_UMAP$UMAP_2,
       pch=20,col='#C59739',cex=0.3)
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

png('PtrLch_trying_colors.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(projection_UMAP$UMAP_1,
     projection_UMAP$UMAP_2,
     pch=20,cex=0.5,
     col=with(color_table,
              Color[match(projection_UMAP$Seurat_clusters,Seurat_cluster)]),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrLch')
legend(10.5,9,xpd=T,bty='n',
       pch=20,col=color_table$Color,
       legend=color_table$Seurat_cluster,pt.cex=2)
for(i in 1:nrow(color_table)) text(Seurat_cluster_center[1,i],
                                   Seurat_cluster_center[2,i],
                                   color_table$Seurat_cluster[i],cex=1.5)
dev.off()

##Change some color by hand
handcraft_color_table = color_table
handcraft_color_table[16,'Color'] = '#2AA02A'
handcraft_color_table[8,'Color'] = '#FF7F0F'
handcraft_color_table[6,'Color'] = '#951DA2'
handcraft_color_table[12,'Color'] = '#FF3E3F'
handcraft_color_table[3,'Color'] = '#8C6447'
handcraft_color_table[c(23,32),'Color'] = '#9B9B9B'

png('PtrLch_Both_Seurat_UMAP_colored_by_Seurat_cluster.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(projection_UMAP$UMAP_1,
     projection_UMAP$UMAP_2,
     pch=20,cex=0.3,
     col=with(handcraft_color_table,
              Color[match(projection_UMAP$Seurat_clusters,Seurat_cluster)]),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrLch')
# for(i in 1:nrow(color_table)) text(Seurat_cluster_center[1,i],
#                                    Seurat_cluster_center[2,i],
#                                    color_table$Seurat_cluster[i],cex=1.5)
dev.off()

##Color with the Cre-picking color
Cre30_color_table = data.frame(Seurat_cluster = 1:33,
                               Color = c('#BDC2C5','#DF35DB','#FFFF55','#FE3BAC','#D6DCE5',
                                         '#4F3098','#AEE9F4','#4F9BF9','#148F7C','#DF2036',
                                         '#44BF77','#C9DD31','#5C4C32','#2695EA','#FFE21C',
                                         '#99CC00','#69D8E9','#87ABFD','#FE4496','#6DC6FD',
                                         '#F0AB2C','#ADEB4D','#FE7D1C','#16EE91','#FD3A5C',
                                         '#A74831','#951DA2','#BB986A','#282DE8','#595959',
                                         '#FF71AD','#FF3E3F','#A262F0'))
png('PtrLch_Both_Seurat_UMAP_colored_by_Seurat_cluster_Cre30.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(projection_UMAP$UMAP_1,
     projection_UMAP$UMAP_2,
     pch=20,cex=0.3,
     col=with(Cre30_color_table,
              Color[match(projection_UMAP$Seurat_clusters,Seurat_cluster)]),
     xlab='UMAP_1',ylab='UMAP_2',main='PtrLch')
for(i in 1:nrow(color_table)) text(Seurat_cluster_center[1,i],
                                   Seurat_cluster_center[2,i],
                                   color_table$Seurat_cluster[i],cex=2.5)
dev.off()


#Conduct the differential expression analysis (DEA) for each handcraft-color-cluster
projection_UMAP$Seurat_color = with(handcraft_color_table,
                                    Color[match(as.numeric(Combined_object@meta.data$seurat_clusters),
                                                Seurat_cluster)])
saveRDS(projection_UMAP,'RDS_PtrLch_Seurat_projection_UMAP.rds')
dim(projection_UMAP) # 7682 8

Combined_object = SetIdent(Combined_object,
                           value=projection_UMAP$Seurat_color)
##One-versus-all DEA
all.markers = FindAllMarkers(object=Combined_object, test.use='MAST', only.pos=T)
all.markers %>% str # 4392 7
saveRDS(all.markers,'RDS_PtrLch_allmarkers.rds')

##Take two for checks
filter(all.markers,cluster=='#D62728') %>% head
FeaturePlot(Combined_object,features='Cluster-1885',pt.size=1)
FeaturePlot(Combined_object,features='Cluster-435',pt.size=1)



# library(GSAR)
# projection_PCA = Combined_object@reductions$pca@cell.embeddings %>% t
# group_ID = ifelse(substr(colnames(projection_PCA),1,3)=='Ptr',1,2)
# result_WW = WWtest(projection_PCA, group_ID, nperm=1000, pvalue.only=F)
# saveRDS(result_WW,'RDS_PtrLch_WWtestResult.rds')



# (鑑定fiber)
# 39
# 210
# 294
# 757
# 988
# 3710
# 6098
# (鑑定organizer)
# 340 
# (鑑定fusiform)
# 487
# 698
# 2653
# 2783
# 4004
# 11791

selected_barcodes = Lch_projection_UMAP$Barcode
Lch_UMI_matrix = readRDS('RDS_Lch_SCseq_UMI_matrix.rds')[,selected_barcodes]
Lch_PtrLch_UMI_matrix = readRDS('RDS_Lch_SCseq_ortho_cluster_UMI_matrix.rds')[,selected_barcodes]

Lch_norm_factor = 1000/apply(Lch_UMI_matrix,2,sum)
Lch_PtrLch_norm_UMI_matrix = apply(Lch_PtrLch_UMI_matrix,1,multiply_by,Lch_norm_factor) %>% t

for(CL in c(39,210,294,757,988,3710,6098,340,487,698,2653,2783,4004,11791)){
    png(paste0('PtrLch_Lch_Seurat_UMAP_colored_with_Cluster',CL,'.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot_log2_UMI = log2(Lch_PtrLch_norm_UMI_matrix[paste0('Cluster_',CL),]+1)
    col_order = order(plot_log2_UMI)
    plot_col_ind = round(plot_log2_UMI/max(plot_log2_UMI)*99)+1
    COL = c('#c8c8c8',colorRampPalette(c('#fdd49e','#ed6044','#7f0000'))(99))
    plot(0,0,type='n',
         xlab='',ylab='',
         xlim=range(Lch_projection_UMAP$UMAP_1),
         ylim=range(Lch_projection_UMAP$UMAP_2),axes=F)
    points(Lch_projection_UMAP$UMAP_1[col_order],
           Lch_projection_UMAP$UMAP_2[col_order],
           col=COL[plot_col_ind][col_order],pch=20,cex=0.5)
    # ylim = c(2.5,5.5)
    # for(i in seq_along(COL)){
    #     x = rep(-8,2)
    #     y = seq(ylim[1],ylim[2],
    #             length.out=length(COL)+1)[c(i,i+1)]
    #     lines(x,y,lwd=20,lend='butt',col=COL[i],xpd=T)
    # }
    # text(-7.5,ylim[2]-0.1,sprintf('%.1f',max(plot_log2_UMI)),adj=0)
    # text(-7.5,ylim[1]+0.1,sprintf('%.1f',min(plot_log2_UMI)),adj=0)
    # text(-8.4,ylim[2]+0.3,'Log2 Exp',adj=0)
    dev.off()
}

png('PtrLch_Lch_Seurat_UMAP_colored_by_CellRanger_cluster_withoutFrame.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot(0,0,type='n',las=1,
     xlim=range(Lch_projection_UMAP$UMAP_1),
     ylim=range(Lch_projection_UMAP$UMAP_2),
     xlab='',ylab='',axes=F)
points(Lch_projection_UMAP$UMAP_1,
       Lch_projection_UMAP$UMAP_2,
       pch=20,cex=0.5,
       col=Lch_projection_UMAP$CellRanger_color)
dev.off()




Ptr = SetIdent(Ptr,value=with(Ptr_barcode2color,
                              Color[match(sub('Ptr_','',rownames(Ptr@meta.data)),Barcode)]))
Ptr_all.markers = FindAllMarkers(object=Ptr, test.use='MAST', only.pos=T)
saveRDS(Ptr_all.markers,'RDS_PtrLch_Ptr_allmarkers.rds')

Lch = SetIdent(Lch,value=with(Lch_barcode2color,
                              Color[match(sub('Lch_','',rownames(Lch@meta.data)),Barcode)]))
Lch_all.markers = FindAllMarkers(object=Lch, test.use='MAST', only.pos=T)
saveRDS(Lch_all.markers,'RDS_PtrLch_Lch_allmarkers.rds')

