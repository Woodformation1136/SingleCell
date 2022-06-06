library(dplyr)
library(slingshot)
library(magrittr)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210609_4S_Pseudoline')

Seurat_projection_UMAP = readRDS('RDS_4S_Seurat_projection_UMAP.rds')


plot_lineage_on_Seurat_UMAP = function(SP='Ptr'){
    
    plot_Seurat_projection_UMAP = filter(Seurat_projection_UMAP,Species==SP)
    
    if(SP=='Lch'){
        foo = plot_Seurat_projection_UMAP %>%
            filter(CellRanger_color=='#9467BD') %>% arrange(UMAP_1)
        excluded_barcodes = foo$SpBarcode %>% head(floor(nrow(foo)*5/100))
        plot_Seurat_projection_UMAP %<>% filter(!(SpBarcode%in%excluded_barcodes))
    }
    
    #Identify the psuedotime of each cell on each lineage
    lineage = new('SlingshotDataSet')
    # str(lineage)
    
    ##Prepare the reducedDim
    rd = as.matrix(plot_Seurat_projection_UMAP[,1:2])
    
    ##Prepare the clusterLabels
    cluster_factor = factor(plot_Seurat_projection_UMAP$CellRanger_color)
    # plot(1:10,1:10,col=levels(cluster_factor),pch=20,cex=5)
    cl = sapply(as.numeric(cluster_factor),
                function(i){
                    out = rep(0,10); out[i] = 1
                    return(out)
                }) %>% t
    rownames(cl) = rownames(rd)
    colnames(cl) = 1:10
    
    ##Prepare the lineages
    lin = sapply(list(c('#FF7F0F','#F8E71C','#E377C2'),
                      c('#9467BD','#2AA02A','#8C564B','#1F77B4'),
                      c('#9467BD','#2AA02A','#8C564B','#D62728')),
                 function(L)
                     na.omit(as.character(match(L,levels(cluster_factor))))) %>%
        set_names(c('Lineage1','Lineage2','Lineage3'))
    
    ##Prepare the adjacency
    adc = matrix(0,10,10)
    rownames(adc) = 1:10
    colnames(adc) = 1:10
    for(L in lin){
        for(Ci in seq(length(L)-1)){
            adc[L[Ci],L[Ci+1]] = 1
            adc[L[Ci+1],L[Ci]] = 1
        }
    }
    
    ##Fill in the contents
    lineage@reducedDim = rd
    lineage@clusterLabels = cl
    lineage@lineages = lin
    lineage@adjacency = adc
    
    ##Plot the lineages
    # plot(plot_Seurat_projection_UMAP$UMAP_1,
    #      plot_Seurat_projection_UMAP$UMAP_2,
    #      pch=20,cex=0.3,
    #      col=plot_Seurat_projection_UMAP$CellRanger_color,
    #      xlab='UMAP_1',ylab='UMAP_2',main='PtrEgr')
    # lines(lineage, lwd=3, col='black')
    
    ##Get lineage curves
    lineage = getCurves(lineage,extend='n')
    for(i in 1:3){
        png(paste0(SP,'_Seurat_UMAP_with_lineage',i,'.png'),
            pointsize=10,width=20,height=15,units='cm',res=300)
        plot(plot_Seurat_projection_UMAP$UMAP_1,
             plot_Seurat_projection_UMAP$UMAP_2,
             pch=20,cex=0.3,
             col=plot_Seurat_projection_UMAP$CellRanger_color,
             xlab='UMAP_1',ylab='UMAP_2',main=SP)
        lines(lineage@curves[[i]], lwd=3, col='black')
        dev.off()
    }
    png(paste0(SP,'_Seurat_UMAP_with_lineage.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    plot(plot_Seurat_projection_UMAP$UMAP_1,
         plot_Seurat_projection_UMAP$UMAP_2,
         pch=20,cex=0.3,
         col=plot_Seurat_projection_UMAP$CellRanger_color,
         xlab='UMAP_1',ylab='UMAP_2',main=SP)
    lines(lineage, lwd=3, col='black')
    dev.off()
}

plot_lineage_on_Seurat_UMAP('Ptr')
plot_lineage_on_Seurat_UMAP('Egr')
plot_lineage_on_Seurat_UMAP('Tar')
plot_lineage_on_Seurat_UMAP('Lch')



