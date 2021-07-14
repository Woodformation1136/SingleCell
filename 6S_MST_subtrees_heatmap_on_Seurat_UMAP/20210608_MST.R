library(igraph)
library(magrittr)
library(Matrix)
library(MASS)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210608_MST_subtrees_heatmap')

PtrEgr_combinedObject = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Egr.rds')
PtrTar_combinedObject = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Tar.rds')
PtrLch_combinedObject = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Lch.rds')
PtrAth2019R_combinedObject = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2019R.rds')
PtrAth2021L_combinedObject = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2021L.rds')
PtrOsa2021R_combinedObject = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Osa2021R.rds')

# PtrEgr_projectionPCA = PtrEgr_combinedObject@reductions$pca@cell.embeddings
# PtrTar_projectionPCA = PtrTar_combinedObject@reductions$pca@cell.embeddings
# PtrLch_projectionPCA = PtrLch_combinedObject@reductions$pca@cell.embeddings
# PtrAth2019R_projectionPCA = PtrAth2019R_combinedObject@reductions$pca@cell.embeddings
# PtrAth2021L_projectionPCA = PtrAth2021L_combinedObject@reductions$pca@cell.embeddings
# PtrOsa2021R_projectionPCA = PtrOsa2021R_combinedObject@reductions$pca@cell.embeddings

PtrEgr_projectionUMAP = PtrEgr_combinedObject@reductions$umap@cell.embeddings
PtrTar_projectionUMAP = PtrTar_combinedObject@reductions$umap@cell.embeddings
PtrLch_projectionUMAP = PtrLch_combinedObject@reductions$umap@cell.embeddings
PtrAth2019R_projectionUMAP = PtrAth2019R_combinedObject@reductions$umap@cell.embeddings
PtrAth2021L_projectionUMAP = PtrAth2021L_combinedObject@reductions$umap@cell.embeddings
PtrOsa2021R_projectionUMAP = PtrOsa2021R_combinedObject@reductions$umap@cell.embeddings

# getMSTsubtreeCenter = function(projection){
#     print('Calculate the distance between each pair of cells')
#     distMatrix = as.matrix(dist(projection)) %>% Matrix(sparse=T)
#     stopifnot(sum(distMatrix==0) == nrow(projection))
# 
#     print('Create the graph from adjacent matrix')
#     graphFull = graph_from_adjacency_matrix(distMatrix,
#                                             mode='undirected',weighted=T)
# 
#     print('Construct the MST')
#     graphMST = mst(graphFull)
# 
#     print('Remove inter-species edges')
#     edgeVname = attr(E(graphMST),'vnames')
#     delEdge = edgeVname %>% grep('Ptr',.,value=T) %>% grep('Egr|Tar|Lch|Ath|Osa',.,value=T)
#     graphCutMST = delete_edges(graphMST,delEdge)
# 
#     print('Extract the subgraph centers')
#     subgraphCenter = c()
#     candidateVertices = attr(V(graphCutMST),'name')
#     while(length(candidateVertices)>0){
#         pickedVertex = candidateVertices[1]
#         pickedVertices = attr(subcomponent(graphCutMST,pickedVertex),'name')
#         pickedGraph = induced_subgraph(graphCutMST,pickedVertices)
#         pickedCloseness = closeness(pickedGraph)
#         if(length(pickedCloseness)==1){
#             subgraphCenter %<>% c(names(pickedCloseness))
#         }else{
#             subgraphCenter %<>% c(names(which.max(pickedCloseness)))
#         }
#         candidateVertices %<>% setdiff(pickedVertices)
#     }
#     return(subgraphCenter)
# }
# 
# print('PtrEgr')
# PtrEgr_subtreeCenter = getMSTsubtreeCenter(PtrEgr_projectionUMAP)
# saveRDS(PtrEgr_subtreeCenter,'RDS_PtrEgr_subtreeCenter.rds')
# 
# print('PtrTar')
# PtrTar_subtreeCenter = getMSTsubtreeCenter(PtrTar_projectionUMAP)
# saveRDS(PtrTar_subtreeCenter,'RDS_PtrTar_subtreeCenter.rds')
# 
# print('PtrLch')
# PtrLch_subtreeCenter = getMSTsubtreeCenter(PtrLch_projectionUMAP)
# saveRDS(PtrLch_subtreeCenter,'RDS_PtrLch_subtreeCenter.rds')
# 
# print('PtrAth2019R')
# PtrAth2019R_subtreeCenter = getMSTsubtreeCenter(PtrAth2019R_projectionUMAP)
# saveRDS(PtrAth2019R_subtreeCenter,'RDS_PtrAth2019R_subtreeCenter.rds')
# 
# print('PtrAth2021L')
# PtrAth2021L_subtreeCenter = getMSTsubtreeCenter(PtrAth2021L_projectionUMAP)
# saveRDS(PtrAth2021L_subtreeCenter,'RDS_PtrAth2021L_subtreeCenter.rds')
# 
# print('PtrOsa2021R')
# Ptr_ID = PtrOsa2021R_projectionPCA %>% rownames %>% grep('Ptr',.,value=T)
# Osa_ID = PtrOsa2021R_projectionPCA %>% rownames %>% grep('Osa',.,value=T) %>% sample(20000)
# PtrOsa2021R_subtreeCenter = getMSTsubtreeCenter(PtrOsa2021R_projectionPCA[c(Ptr_ID,Osa_ID),])
# saveRDS(PtrOsa2021R_subtreeCenter,'RDS_PtrOsa2021R_subtreeCenter.rds')


PtrEgr_subtreeCenter = readRDS('RDS_PtrEgr_subtreeCenter.rds')
PtrTar_subtreeCenter = readRDS('RDS_PtrTar_subtreeCenter.rds')
PtrLch_subtreeCenter = readRDS('RDS_PtrLch_subtreeCenter.rds')
PtrAth2019R_subtreeCenter = readRDS('RDS_PtrAth2019R_subtreeCenter.rds')
PtrAth2021L_subtreeCenter = readRDS('RDS_PtrAth2021L_subtreeCenter.rds')
PtrOsa2021R_subtreeCenter = readRDS('RDS_PtrOsa2021R_subtreeCenter.rds')


centerContourPlot = function(subgraphCenter = PtrEgr_subtreeCenter,
                             projectionUMAP = PtrEgr_projectionUMAP,
                             lowerPercentile = 20,
                             higherPercentile = 50,
                             main = 'PtrEgr',
                             plot_name = 'PtrEgr'){
    
    Ptr_projectionUMAP = projectionUMAP %>% extract(grepl('Ptr',rownames(.)),)
    Ptr_density_map = kde2d(Ptr_projectionUMAP[,1],
                            Ptr_projectionUMAP[,2],n=500,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
    get_territory_density = function(Coor){
        out = Ptr_density_map$z[max(which(Ptr_density_map$x < Coor[1])),
                                max(which(Ptr_density_map$y < Coor[2]))]
        return(out)
    }
    Ptr_territory_density = Ptr_projectionUMAP %>% apply(1,get_territory_density)
    norm_factor = 1/sum(Ptr_territory_density)
    # print(norm_factor)
    
    Ptr_SubgraphCenter = subgraphCenter %>% extract(grepl('Ptr',.))
    Sp2_SubgraphCenter = subgraphCenter %>% extract(!grepl('Ptr',.))
    
    num_Ptr =  sum(grepl('Ptr',rownames(projectionUMAP)))
    num_Sp2 =  sum(!grepl('Ptr',rownames(projectionUMAP)))
    num_total = nrow(projectionUMAP)
    num_center_Ptr = length(Ptr_SubgraphCenter)
    num_center_Sp2 = length(Sp2_SubgraphCenter)

    get_subgraphcenter_density = function(partSubgraphCenter){
        allCenterCoordinate = projectionUMAP[subgraphCenter,]
        partCenterCoordinate = projectionUMAP[partSubgraphCenter,]
        density_map = kde2d(partCenterCoordinate[,1],
                            partCenterCoordinate[,2],n=500,
                            h=apply(allCenterCoordinate,2,bandwidth.nrd)/2,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
        return(density_map)
    }
    Ptr_subgraphcenter_density = get_subgraphcenter_density(Ptr_SubgraphCenter)
    Sp2_subgraphcenter_density = get_subgraphcenter_density(Sp2_SubgraphCenter)
    stopifnot(Ptr_subgraphcenter_density$x==Sp2_subgraphcenter_density$x)
    stopifnot(Ptr_subgraphcenter_density$y==Sp2_subgraphcenter_density$y)
    
    merge_subgraphcenter_density = Ptr_subgraphcenter_density
    merge_subgraphcenter_density$z =
        (num_Ptr/num_total)*num_center_Ptr*Ptr_subgraphcenter_density$z +
        (num_Sp2/num_total)*num_center_Sp2*Sp2_subgraphcenter_density$z
    
    merge_subgraphcenter_density$z %<>% multiply_by(norm_factor)
    print(paste0('Plot density max:',max(merge_subgraphcenter_density$z)))
    print(paste0('Plot density Q75:',quantile(merge_subgraphcenter_density$z,0.75)))
    print(paste0('Plot density Q50:',quantile(merge_subgraphcenter_density$z,0.50)))
    print(paste0('Plot density Q25:',quantile(merge_subgraphcenter_density$z,0.25)))
    
    territory_map = kde2d(projectionUMAP[,1],
                          projectionUMAP[,2],n=500,h=0.02,
                          lims=c(min(projectionUMAP[,1])-0.5,
                                 max(projectionUMAP[,1])+0.5,
                                 min(projectionUMAP[,2])-0.5,
                                 max(projectionUMAP[,2])+0.5))
    
    png(paste0(plot_name,'.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    {
        plot(NA,
             xlim=range(projectionUMAP[,'UMAP_1']),
             ylim=range(projectionUMAP[,'UMAP_2']),
             xlab='',ylab='',axes=F,main=main)
        # Levels = seq(0,max(density_map$z),length.out=200)
        Levels = seq(0,1,length.out=500)
        lowerRank = 500 * lowerPercentile/100
        higherRank = 500 * higherPercentile/100
        .filled.contour(merge_subgraphcenter_density$x,
                        merge_subgraphcenter_density$y,
                        merge_subgraphcenter_density$z,
                        levels=Levels,
                        col=c(colorRampPalette(c('#58b94e','#fec84a'))(lowerRank),
                              colorRampPalette(c('#fec84a','#eb473e'))(higherRank-lowerRank),
                              colorRampPalette(c('#eb473e','#713020'))(500-higherRank)))
        # col=c(colorRampPalette(c('#dde6ff','#ffffff'))(whiteRank),
        #       colorRampPalette(c('#ffffff','#ff0000'))(redRank-whiteRank),
        #       heat.colors(500-redRank))
        .filled.contour(territory_map$x,
                        territory_map$y,
                        ifelse(territory_map$z>0,1,0),
                        levels=c(0,0.5),col=c('white',NA))
        contour(territory_map$x,territory_map$y,ifelse(territory_map$z>0,1,0),
                levels=0.5,lwd=0.8,drawlabels=F,add=T)
        # COL = rep('gray90',nrow(projectionUMAP))
        # COL[match(subgraphCenter,rownames(projectionUMAP))] = 2
        # colInd = order(COL,decreasing=T)
        # points(projectionUMAP[colInd,'UMAP_1'],
        #        projectionUMAP[colInd,'UMAP_2'],
        #        col=COL[colInd],pch=20,cex=1)
    }
    dev.off()
}


SPs = c('Egr','Tar','Lch','Ath2019R','Ath2021L','Osa2021R')
for(SP in SPs){
    print(SP)
    code = paste0("subgraphCenter = Ptr",SP,"_subtreeCenter")
    eval(parse(text=code))
    code = paste0("projectionUMAP = Ptr",SP,"_projectionUMAP")
    eval(parse(text=code))
    
    centerContourPlot(subgraphCenter,
                      projectionUMAP,
                      lowerPercentile = 5,
                      higherPercentile = 40,
                      main = paste0('Ptr',SP),
                      plot_name = paste0('Ptr',SP))
}
