setwd('C:\\Users\\Jung-Chen\\Desktop\\Single_cell\\20210616_Ptr_CellRanger_UMAP')

Ptr_SCseq_cell_cluster_10 = readRDS('RDS_Ptr_SCseq_cell_cluster_10.rds')
Ptr_projection_umap_CellRanger = readRDS('RDS_Ptr_projection_umap.rds')

COL = c('#1F77B4','#8C564B','#FF7F0F','#2AA02A','#F8E71C',
        '#9467BD','#D62728','#E377C2','#9B9B9B','#4B4B4B')
png('Ptr_CellRanger_UMAP.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
par(mar=c(0,0,0,0))
plot(Ptr_projection_umap_CellRanger$UMAP.1,
     Ptr_projection_umap_CellRanger$UMAP.2,
     pch=20,cex=1.2,
     col=COL[Ptr_SCseq_cell_cluster_10$Cluster])
dev.off()