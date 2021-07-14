library(Seurat)
library(magrittr)
library(filesstrings)
library(scales)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210526_10x_format_Ptr_Ath\\2019_Ryu')


#Prepare gz files for input
all_files = list.files(recursive=T)
all_dirs = list.dirs()

bar_vector = grep('barcodes.tsv$',all_files,value=T)
for(i in seq_along(bar_vector)){
    bar = bar_vector[i]
    out = gzfile(paste0(bar,'.gz'),'w')
    cat(paste0('Ath_2019R',i,'_',scan(bar,what='')),
        sep='\n',file=out)
    close(out)
    file.move(paste0(bar,'.gz'),
              sub('barcodes.tsv.gz','processed',paste0(bar,'.gz')),
              overwrite=T)
}
for(bar in grep('genes.tsv$',all_files,value=T)){
    system(paste('gzip --force --keep',bar))
    file_path = sub('genes.tsv','',bar)
    file.move(paste0(file_path,'genes.tsv.gz'),
              paste0(file_path,'processed'),
              overwrite=T)
    file.rename(paste0(file_path,'processed/genes.tsv.gz'),
                paste0(file_path,'processed/features.tsv.gz'))
}
for(bar in grep('matrix.mtx$',all_files,value=T)){
    system(paste('gzip --force --keep',bar))
    file_path = sub('matrix.mtx','',bar)
    file.move(paste0(file_path,'matrix.mtx.gz'),
              paste0(file_path,'processed'),
              overwrite=T)
}


#Setup the Seurat Object
Ath2019R1.data = Read10X(data.dir = 'Sample_WT-WERGFP\\filtered_gene_bc_matrices\\TAIR10\\processed')
Ath2019R2.data = Read10X(data.dir = 'Sample_WT-WERGFP_2\\filtered_gene_bc_matrices\\TAIR10\\processed')
Ath2019R3.data = Read10X(data.dir = 'Sample_WT-WERGFP_3\\filtered_gene_bc_matrices\\TAIR10\\processed')
Ath2019R1 = CreateSeuratObject(counts = Ath2019R1.data, project = 'Ath2019R1', min.cells = 3, min.features = 200)
Ath2019R2 = CreateSeuratObject(counts = Ath2019R2.data, project = 'Ath2019R2', min.cells = 3, min.features = 200)
Ath2019R3 = CreateSeuratObject(counts = Ath2019R3.data, project = 'Ath2019R3', min.cells = 3, min.features = 200)
Ath2019R1
Ath2019R2
Ath2019R3

#Normalize the data
Ath2019R1 = NormalizeData(Ath2019R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R2 = NormalizeData(Ath2019R2, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R3 = NormalizeData(Ath2019R3, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ath2019R1 = FindVariableFeatures(Ath2019R1, selection.method = 'vst', nfeatures = 2000)
Ath2019R2 = FindVariableFeatures(Ath2019R2, selection.method = 'vst', nfeatures = 2000)
Ath2019R3 = FindVariableFeatures(Ath2019R3, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Ath2019R = FindIntegrationAnchors(object.list = list(Ath2019R1,Ath2019R2,Ath2019R3),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Ath2019R = IntegrateData(anchorset = integration_anchors_Ath2019R,
                                         preserve.order = T,
                                         sample.tree = rbind(c(-1,-2),
                                                             c(1,-3)))


#Run the standard workflow for visualization and clustering
Combined_object_Ath2019R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Ath2019R %<>% FindNeighbors(reduction='pca', dims=1:30) %>% FindClusters(resolution=0.5)

saveRDS(Combined_object_Ath2019R,'RDS_Combined_object_SF10000_AF2000_KA5_Ath2019R.rds')
Combined_object_Ath2019R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ath2019R.rds')

# DimPlot(Combined_object_Ath2019R,
#         reduction = 'umap',
#         group.by = 'orig.ident',
#         pt.size = 0.5)
# DimPlot(Combined_object_Ath2019R,
#         reduction = 'umap',
#         group.by = 'seurat_clusters',
#         pt.size = 0.5)
# FeaturePlot(Combined_object_Ath2019R,
#             features = 'BHLH30',
#             pt.size = 1, order = T,
#             slot = 'counts')


projection_UMAP = Combined_object_Ath2019R@reductions$umap@cell.embeddings %>% as.data.frame



png('batch.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot_batch_id = Combined_object_Ath2019R@meta.data %>% 
    rownames %>% strsplit('_') %>% sapply(extract,2) %>% sub('2019R','',.) %>% as.numeric
COL = hue_pal()(3)
plot(0,0,type='n',
     xlab='',ylab='',
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),axes=F)
points(projection_UMAP$UMAP_1,
       projection_UMAP$UMAP_2,
       col=COL[plot_batch_id],pch=20,cex=0.5)
legend(-15,20,bty='n',
       legend=paste('Replicate',1:3),col=COL,
       pch=15,pt.cex=2,y.intersp=1.2,xpd=T)
dev.off()



# foo = read.table('Sample_WT-WERGFP/filtered_gene_bc_matrices/TAIR10/genes.tsv',
#                  header=F,sep='\t')
# with(foo,V2[V1=='AT1G68810']) # 'BHLH30'

png('1.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot_log2_UMI = log2(Combined_object_Ath2019R@assays$RNA@counts['BHLH30',]+1)
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
ylim = c(5,15)
for(i in seq_along(COL)){
    x = rep(-12,2)
    y = seq(ylim[1],ylim[2],
            length.out=length(COL)+1)[c(i,i+1)]
    lines(x,y,lwd=20,lend='butt',col=COL[i],xpd=T)
}
text(-11,ylim[2]-0.4,sprintf('%.1f',max(plot_log2_UMI)),adj=0)
text(-11,ylim[1]+0.4,sprintf('%.1f',min(plot_log2_UMI)),adj=0)
text(-12.5,ylim[2]+2,'Log2 Exp\nBHLH30',adj=0)
dev.off()



png('c.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
plot_cluster_id = as.numeric(Combined_object_Ath2019R@meta.data$seurat_clusters)
COL = hue_pal()(16); COL[10] = '#d0021b'
plot(0,0,type='n',
     xlab='',ylab='',
     xlim=range(projection_UMAP$UMAP_1),
     ylim=range(projection_UMAP$UMAP_2),axes=F)
points(projection_UMAP$UMAP_1,
       projection_UMAP$UMAP_2,
       col=COL[plot_cluster_id],pch=20,cex=0.5)
legend(-15,20,bty='n',
       legend=paste('Cluster',1:16),col=COL,
       pch=15,pt.cex=2,y.intersp=1.2,xpd=T)
dev.off()


selected_barcodes = rownames(Combined_object_Ath2019R@meta.data)[plot_cluster_id==10]
length(selected_barcodes) # 353
cat(selected_barcodes,file='Xylem_barcodes.txt',sep='\n')


