library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
oriPar = par(no.readonly=T)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210607_Dotplot_with_SCseq_DEG')

Ptr_SCseq_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')
Ptr_SCseq_norm_UMI_matrix = apply(Ptr_SCseq_UMI_matrix,2,function(x)1000*x/sum(x))
Ptr_SCseq_cell_cluster_10 = readRDS('RDS_Ptr_SCseq_cell_cluster_10.rds')

Ptr_DEG  = list()
for(i in 1:8){
    code = paste0("Ptr_DEG_table = read.csv('Ptr0",i,".csv')")
    eval(parse(text=code))
    
    Ptr_DEG[[i]] = Ptr_DEG_table %>%
        filter(Mean.UMI.count>0.5,Log2.fold.change>2,Adjusted.p.value<0.05) %>%
        arrange(-Log2.fold.change) %>% extract(,'Gene.ID') %>% as.character()
}



marker_dotplot = function(plot_geneID,
                          plot_name = 'Dotplot.png',width=15,height=15){
    scale_mean_exp_long = sapply(1:10,
                                 function(cl){
                                     selected_barcode = with(Ptr_SCseq_cell_cluster_10,
                                                             Barcode[Cluster==cl])
                                     mean_exp = Ptr_SCseq_norm_UMI_matrix %>%
                                         extract(plot_geneID,selected_barcode,drop=F) %>%
                                         apply(1,mean)
                                     return(mean_exp)
                                 }) %>% 
        apply(1,function(x)x/max(x)*100) %>% t %>% 
        as.data.frame %>% 
        set_colnames(paste0('cellCluster_',1:10)) %>%
        inset(,'geneID',rownames(.)) %>%
        melt(id.vars='geneID',
             variable.name='cellCluster',
             value.name='scaleMeanExp')
    
    prop_exp_long = sapply(1:10,
                           function(cl){
                               selected_barcode = with(Ptr_SCseq_cell_cluster_10,
                                                       Barcode[Cluster==cl])
                               prop_exp = Ptr_SCseq_norm_UMI_matrix %>%
                                   extract(plot_geneID,selected_barcode,drop=F) %>%
                                   apply(1,function(x)100*sum(x>1)/length(x))
                               return(prop_exp)
                           }) %>% 
        as.data.frame %>%
        set_colnames(paste0('cellCluster_',1:10)) %>%
        inset(,'geneID',rownames(.)) %>%
        melt(id.vars='geneID',
             variable.name='cellCluster',
             value.name='propExp')
    
    stopifnot(all(scale_mean_exp_long$geneID == prop_exp_long$geneID))
    stopifnot(all(scale_mean_exp_long$cellCluster == prop_exp_long$cellCluster))
    
    plot_long_table = scale_mean_exp_long
    plot_long_table$propExp = prop_exp_long$propExp
    plot_long_table$geneID %<>% sub('Potri.','',.) %>% sub('.v4.1','',.)
    # dim(plot_long_table) # all: 160 4
    # head(plot_long_table)
    
    p1 = ggplot(plot_long_table) +
        geom_point(aes(x=cellCluster,
                       y=factor(geneID,
                                levels=plot_geneID%>%sub('Potri.','',.)%>%sub('.v4.1','',.)%>%rev),
                       size=propExp,
                       color=scaleMeanExp)) +
        xlab('Cell cluster') + ylab('Gene ID') +
        scale_x_discrete(breaks=paste0('cellCluster_',1:10),
                         labels=1:10) +
        scale_size_area(limits=c(0,100),breaks=seq(25,100,25)) +
        scale_color_gradient(low='gray90',high='black',na.value='gray90',
                             limit=c(0,100),breaks=seq(0,100,50)) +
        labs(size='Cell expressed (%)',col='Expression') +
        theme_classic() +
        theme(axis.text.x=element_text(angle=0, vjust=0.5),
              legend.key.size = unit(0.03,'npc'),
              legend.position='right')
    
    ggsave(plot_name,
           plot=p1,
           device='png',
           width=width,height=height,
           units='cm',dpi=300)
}


# for(i in 1:8){
#     plot_geneIDs = na.omit(Ptr_DEG[[i]][1:100])
#     marker_dotplot(plot_geneIDs,
#                    paste0('Dotplot_',i,'.png'),
#                    width=16,height=60*(length(plot_geneIDs)/100))
# }


# Cluster1
C1 = c('Potri.005G233600.v4.1',
       'Potri.008G020900.v4.1',
       'Potri.017G015700.v4.1',
       'Potri.019G069300.v4.1')

# Cluster2
C2 = c('Potri.002G218725.v4.1',
       'Potri.005G081400.v4.1',
       'Potri.006G053300.v4.1',
       'Potri.014G022200.v4.1')

# Cluster3
C3 = c('Potri.005G166100.v4.1',
       'Potri.008G043900.v4.1')

# Cluster4
C4 = c('Potri.005G007200.v4.1',
       'Potri.007G020900.v4.1',
       'Potri.008G225001.v4.1')

# Cluster5
C5 = c('Potri.008G174100.v4.1',
       'Potri.018G114300.v4.1')

# Cluster6
C6 = c('Potri.006G196900.v4.1',
       'Potri.004G002500.v4.1',
       'Potri.001G154100.v4.1',
       'Potri.015G048700.v4.1')

# Cluster7
C7 = c('Potri.002G257900.v4.1',
       'Potri.006G129200.v4.1',
       'Potri.010G141600.v4.1',
       'Potri.015G060100.v4.1')

# Cluster8
C8 = c('Potri.002G104600.v4.1',
       'Potri.002G162400.v4.1',
       'Potri.002G251000.v4.1',
       'Potri.006G106900.v4.1')

# 12467
marker_dotplot(c(C1,C2,C4,C6,C7),
               paste0('Dotplot_12467.png'),
               width=16,height=12)
# 358
marker_dotplot(c(C3,C5,C8),
               paste0('Dotplot_358.png'),
               width=16,height=6)
