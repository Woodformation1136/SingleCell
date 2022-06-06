library(Seurat)
library(magrittr)
library(filesstrings)
library(scales)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210527_10x_format_Ptr_Osa\\2021_Zhang')


#Setup the Seurat Object
Osa2021R1.data = Read10X(data.dir = 'osRoot1\\filtered_feature_bc_matrix',gene.column=1)
Osa2021R2.data = Read10X(data.dir = 'osRoot2\\filtered_feature_bc_matrix',gene.column=1)
colnames(Osa2021R1.data) %<>% paste0('Osa_2021R1_',.)
colnames(Osa2021R2.data) %<>% paste0('Osa_2021R2_',.)
Osa2021R1 = CreateSeuratObject(counts = Osa2021R1.data, project = 'Osa2021R1', min.cells = 3, min.features = 200)
Osa2021R2 = CreateSeuratObject(counts = Osa2021R2.data, project = 'Osa2021R2', min.cells = 3, min.features = 200)

#Normalize the data
Osa2021R1 = NormalizeData(Osa2021R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Osa2021R2 = NormalizeData(Osa2021R2, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Osa2021R1 = FindVariableFeatures(Osa2021R1, selection.method = 'vst', nfeatures = 2000)
Osa2021R2 = FindVariableFeatures(Osa2021R2, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Osa2021R = FindIntegrationAnchors(object.list = list(Osa2021R1,Osa2021R2),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Osa2021R = IntegrateData(anchorset = integration_anchors_Osa2021R)


#Run the standard workflow for visualization and clustering
Combined_object_Osa2021R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30, seed.use = 2021)
Combined_object_Osa2021R %<>% FindNeighbors(reduction='pca', dims=1:30, k.param = 10) %>% FindClusters(resolution=0.5)

saveRDS(Combined_object_Osa2021R,'RDS_Combined_object_SF10000_AF2000_KA5_Osa2021R.rds')
Combined_object_Osa2021R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Osa2021R.rds')

# DimPlot(Combined_object_Osa2021R,
#         reduction = 'umap',
#         group.by = 'orig.ident',
#         pt.size = 0.5)
# DimPlot(Combined_object_Osa2021R,
#         reduction = 'umap',
#         group.by = 'seurat_clusters',
#         pt.size = 0.5)
# FeaturePlot(Combined_object_Osa2021R,
#             features = 'Os10g0467800',
#             pt.size = 1, order = T,
#             slot = 'counts')


projection_UMAP = Combined_object_Osa2021R@reductions$umap@cell.embeddings %>% as.data.frame



png('batch.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot_batch_id = Combined_object_Osa2021R@meta.data %>% 
    rownames %>% strsplit('_') %>% sapply(extract,2) %>% sub('2021R','',.) %>% as.numeric
COL = hue_pal()(2)
plot(0,0,type='n',
     xlab='',ylab='',
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),axes=F)
shuffle_id = sample(seq_along(plot_batch_id))
points(projection_UMAP$UMAP_1[shuffle_id],
       projection_UMAP$UMAP_2[shuffle_id],
       col=COL[plot_batch_id][shuffle_id],pch=20,cex=0.5)
legend(-10,-7,bty='n',
       legend=paste('Replicate',1:2),col=COL,
       pch=15,pt.cex=2,y.intersp=1.2,xpd=T)
dev.off()



# foo = read.table('Sample_WT-WERGFP/filtered_gene_bc_matrices/TAIR10/genes.tsv',
#                  header=F,sep='\t')
# with(foo,V2[V1=='AT1G68810']) # 'BHLH30'

xylem_markers = c('Os01g0750300','Os08g0489300','Os07g0638500',
                  'Os10g0467800','Os09g0422500','Os04g0536500')

for(m in 1:6){
    png(paste0(m,'.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot_log2_UMI = log2(Combined_object_Osa2021R@assays$RNA@counts[xylem_markers[m],]+1)
    col_order = order(plot_log2_UMI)
    plot_col_ind = round(plot_log2_UMI/max(plot_log2_UMI)*99)+1
    COL = c('#c8c8c8',colorRampPalette(c('#fdd49e','#ed6044','#7f0000'))(99))
    plot(0,0,type='n',
         xlab='',ylab='',
         xlim=range(projection_UMAP$UMAP_1),
         ylim=range(projection_UMAP$UMAP_2),axes=F)
    points(projection_UMAP$UMAP_1[col_order],
           projection_UMAP$UMAP_2[col_order],
           col=COL[plot_col_ind][col_order],pch=20,cex=0.5)
    ylim = c(0,10)
    for(i in seq_along(COL)){
        x = rep(-13.8,2)
        y = seq(ylim[1],ylim[2],
                length.out=length(COL)+1)[c(i,i+1)]
        lines(x,y,lwd=20,lend='butt',col=COL[i],xpd=T)
    }
    text(-13.3,ylim[2]-0.2,sprintf('%.1f',max(plot_log2_UMI)),adj=0,xpd=T)
    text(-13.3,ylim[1]+0.2,sprintf('%.1f',min(plot_log2_UMI)),adj=0,xpd=T)
    text(-14.3,ylim[2]+1.2,paste0('Log2 Exp\n',xylem_markers[m]),adj=0,xpd=T)
    dev.off()
}




png('c.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot_cluster_id = as.numeric(Combined_object_Osa2021R@meta.data$seurat_clusters)
COL = hue_pal()(23); COL[c(11,20)] = '#d0021b'
plot(0,0,type='n',
     xlab='',ylab='',
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),axes=F)
points(projection_UMAP$UMAP_1,
       projection_UMAP$UMAP_2,
       col=COL[plot_cluster_id],pch=20,cex=0.5)
legend(-14.5,10,bty='n',
       legend=paste('Cluster',1:23),col=COL,
       pch=15,pt.cex=2,y.intersp=1.2,xpd=T)
dev.off()


selected_barcodes = rownames(Combined_object_Osa2021R@meta.data)[plot_cluster_id%in%c(11,20)]
length(selected_barcodes) # 2508
cat(selected_barcodes,file='Xylem_barcodes.txt',sep='\n')


