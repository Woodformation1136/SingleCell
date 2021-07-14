library(magrittr)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210615_MetaCell_with_CellRanger_cluster')

Ptr_MetaCell = readRDS('RDS_sc_df_PtGenes.rds')[,1:3]
dim(Ptr_MetaCell) # 4473 3


##Input the UMAP and clustering data
Ptr_SCseq_cell_cluster_10 = readRDS('RDS_Ptr_SCseq_cell_cluster_10.rds')
Ptr_projection_umap_CellRanger = readRDS('RDS_Ptr_projection_umap.rds')

##Replot the CellRanger UMAP and color with CellRanger clusters and 
##build the barcode2color table
COL = c('#1F77B4','#8C564B','#FF7F0F','#2AA02A','#F8E71C',
        '#9467BD','#D62728','#E377C2','#9B9B9B','#4B4B4B')
plot(Ptr_projection_umap_CellRanger$UMAP.1,
     Ptr_projection_umap_CellRanger$UMAP.2,
     pch=20,cex=1.2,
     col=COL[Ptr_SCseq_cell_cluster_10$Cluster])
Ptr_barcode2color = data.frame(Barcode = Ptr_SCseq_cell_cluster_10$Barcode,
                               Color = COL[Ptr_SCseq_cell_cluster_10$Cluster])


png('Ptr_MetaCell_colored_by_CellRanger_cluster.png',
    pointsize=10,width=20,height=15,units='cm',res=300)
par(mar=c(0,0,0,0))
plot(Ptr_MetaCell$X,Ptr_MetaCell$Y,
     pch=20,
     col=with(Ptr_barcode2color,Color[match(rownames(Ptr_MetaCell),Barcode)]))
par(oriPar)
dev.off()

