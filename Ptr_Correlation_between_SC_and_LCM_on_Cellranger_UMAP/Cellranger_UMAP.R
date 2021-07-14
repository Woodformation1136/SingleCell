library(magrittr)
library(RColorBrewer)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210520_Cellranger_UMAP_Ptr_recolor')

Ptr_projection_umap = readRDS('RDS_Ptr_projection_umap.rds')
Ptr_LCM_mean_TPM_matrix = readRDS('RDS_Ptr_LCM_mean_TPM_matrix.rds')
Ptr_SCseq_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')

geneID = rownames(Ptr_LCM_mean_TPM_matrix)
Ptr_SCseq_UMI_matrix = Ptr_SCseq_UMI_matrix[geneID,]
all(rownames(Ptr_LCM_mean_TPM_matrix) == rownames(Ptr_SCseq_UMI_matrix))

Cre_correlation_matrix = cor(Ptr_SCseq_UMI_matrix,Ptr_LCM_mean_TPM_matrix)
all(rownames(Cre_correlation_matrix) == Ptr_projection_umap$Barcode)


#GSE81077 Ptr LCM data
##Input the label table
Ptr_LCM_label = read.csv('GSE81077/20210326_Inputfilelist.csv')[,c('Condition','Group')]
##Input the TPM of each sample
for(i in seq(nrow(Ptr_LCM_label))){
    inData = read.table(paste0('GSE81077/outputStringtie_geneAbundances/',Ptr_LCM_label$Group[i],'.tsv'),
                        sep='\t',header=T)
    tempData = inData[,'TPM',drop=F]
    rownames(tempData) = inData$Gene.ID
    if(i==1){
        Ptr_LCM_TPM_matrix = tempData
        geneID = rownames(tempData)
    }else{
        Ptr_LCM_TPM_matrix = cbind(Ptr_LCM_TPM_matrix,
                                   tempData[geneID,])
    }
}
colnames(Ptr_LCM_TPM_matrix) = Ptr_LCM_label$Condition
Ptr_GSE81077_LCM_TPM_matrix = as.matrix(Ptr_LCM_TPM_matrix)
saveRDS(Ptr_GSE81077_LCM_TPM_matrix,'RDS_Ptr_GSE81077_LCM_TPM_matrix.rds')

Ptr_GSE81077_LCM_mean_TPM_matrix = apply(Ptr_GSE81077_LCM_TPM_matrix,1,
                                         function(x){
                                             c(fiber = mean(x[1:3]),
                                               tri_celltype = mean(x[4:6]),
                                               vessel = mean(x[7:9]))
                                         }) %>% t

geneID = rownames(Ptr_GSE81077_LCM_mean_TPM_matrix)
Ptr_SCseq_UMI_matrix = Ptr_SCseq_UMI_matrix[geneID,]
all(rownames(Ptr_GSE81077_LCM_mean_TPM_matrix) == rownames(Ptr_SCseq_UMI_matrix))

GSE81077_correlation_matrix = cor(Ptr_SCseq_UMI_matrix,Ptr_GSE81077_LCM_mean_TPM_matrix)
all(rownames(GSE81077_correlation_matrix) == Ptr_projection_umap$Barcode)



for(data_index in 1:2){
    correlation_matrix = list(Cre_correlation_matrix,
                              GSE81077_correlation_matrix)[[data_index]]
    prefix = c('','GSE81077_')[data_index]
    
    for(pp1 in seq(10,70,10)){
        for(pp2 in seq(pp1+10,90,10)){
            COL = colorRampPalette(c(colorRampPalette(c('#DDDDDD','#EEEEEE'))(pp1),
                                     colorRampPalette(c('#EEEEEE','#005c9c'))(pp2-pp1),
                                     rep('#005c9c',100-pp2)))(100)
            #Original blue: #1F77B4
            
            png(paste0(prefix,'Ptr_CellRanger_UMAP_LCM_vessel_',pp1,'_',pp2,'.png'),
                pointsize=10,width=20,height=15,units='cm',res=300)
            par(mar=c(0,0,0,0))
            plot_correlation = correlation_matrix[,'vessel'] %>% as.numeric
            COL_ind = round(plot_correlation/max(plot_correlation)*99+1)
            order_ind = order(plot_correlation)
            plot(Ptr_projection_umap$UMAP.1[order_ind],
                 Ptr_projection_umap$UMAP.2[order_ind],
                 col=COL[COL_ind][order_ind],
                 pch=20,cex=1.4,
                 main='',xlab='UMAP_1',ylab='UMAP_2')
            for(i in seq_along(COL)){
                x = seq(-9,-2,length.out=length(COL)+1)[c(i,i+1)]
                lines(x,rep(-9,2),lwd=40,lend='butt',col=COL[i])
                # y = seq(-10,-5,length.out=length(COL)+1)[c(i,i+1)]
                # lines(rep(-9,2),y,lwd=20,lend='butt',col=COL[i])
                # if(i==pp1) lines(rep(-9.2,2),y,lwd=10,lend='butt')
                # if(i==pp2) lines(rep(-9.2,2),y,lwd=10,lend='butt')
            }
            text(-2,-7.5,sprintf('%.2f',max(plot_correlation)),cex=4)
            text(-9,-7.5,sprintf('%.2f',0),cex=4)
            # text(-9,-4.6,'Cor')
            # text(-8.5,-5,sprintf('%.2f',max(plot_correlation)),adj=0)
            # text(-8.5,-10,sprintf('%.2f',0),adj=0)
            par(mar=c(5,4,4,2)+0.1)
            dev.off()
        }
    }
    
    for(pp in seq(10,90,10)){
        turning_percentile = pp
        COL = colorRampPalette(c(colorRampPalette(c('#DDDDDD','#EEEEEE'))(turning_percentile),
                                 colorRampPalette(c('#EEEEEE','#D62728'))(100-turning_percentile)))(100)
        
        png(paste0(prefix,'Ptr_CellRanger_UMAP_LCM_fiber_',pp,'.png'),
            pointsize=10,width=20,height=15,units='cm',res=300)
        par(mar=c(0,0,0,0))
        plot_correlation = correlation_matrix[,'fiber'] %>% as.numeric
        COL_ind = round(plot_correlation/max(plot_correlation)*99+1)
        order_ind = order(plot_correlation)
        plot(Ptr_projection_umap$UMAP.1[order_ind],
             Ptr_projection_umap$UMAP.2[order_ind],
             col=COL[COL_ind][order_ind],
             pch=20,cex=1.4,
             main='',xlab='UMAP_1',ylab='UMAP_2')
        for(i in seq_along(COL)){
            x = seq(-9,-2,length.out=length(COL)+1)[c(i,i+1)]
            lines(x,rep(-9,2),lwd=40,lend='butt',col=COL[i])
            # y = seq(-10,-5,length.out=length(COL)+1)[c(i,i+1)]
            # lines(rep(-9,2),y,lwd=20,lend='butt',col=COL[i])
            # if(i==pp1) lines(rep(-9.2,2),y,lwd=10,lend='butt')
            # if(i==pp2) lines(rep(-9.2,2),y,lwd=10,lend='butt')
        }
        text(-2,-7.5,sprintf('%.2f',max(plot_correlation)),cex=4)
        text(-9,-7.5,sprintf('%.2f',0),cex=4)
        # text(-9,-4.6,'Cor')
        # text(-8.5,-5,sprintf('%.2f',max(plot_correlation)),adj=0)
        # text(-8.5,-10,sprintf('%.2f',0),adj=0)
        par(mar=c(5,4,4,2)+0.1)
        dev.off()
    }
    
    if(data_index==1){
        for(pp in seq(10,90,10)){
            turning_percentile = pp
            COL = colorRampPalette(c(colorRampPalette(c('#DDDDDD','#EEEEEE'))(turning_percentile),
                                     colorRampPalette(c('#EEEEEE','#E377C2'))(100-turning_percentile)))(100)
            
            png(paste0(prefix,'Ptr_CellRanger_UMAP_LCM_ray_',pp,'.png'),
                pointsize=10,width=20,height=15,units='cm',res=300)
            par(mar=c(0,0,0,0))
            plot_correlation = correlation_matrix[,'ray'] %>% as.numeric
            COL_ind = round(plot_correlation/max(plot_correlation)*99+1)
            order_ind = order(plot_correlation)
            plot(Ptr_projection_umap$UMAP.1[order_ind],
                 Ptr_projection_umap$UMAP.2[order_ind],
                 col=COL[COL_ind][order_ind],
                 pch=20,cex=1.4,
                 main='',xlab='UMAP_1',ylab='UMAP_2')
            for(i in seq_along(COL)){
                x = seq(-9,-2,length.out=length(COL)+1)[c(i,i+1)]
                lines(x,rep(-9,2),lwd=40,lend='butt',col=COL[i])
                # y = seq(-10,-5,length.out=length(COL)+1)[c(i,i+1)]
                # lines(rep(-9,2),y,lwd=20,lend='butt',col=COL[i])
                # if(i==pp1) lines(rep(-9.2,2),y,lwd=10,lend='butt')
                # if(i==pp2) lines(rep(-9.2,2),y,lwd=10,lend='butt')
            }
            text(-2,-7.5,sprintf('%.2f',max(plot_correlation)),cex=4)
            text(-9,-7.5,sprintf('%.2f',0),cex=4)
            # text(-9,-4.6,'Cor')
            # text(-8.5,-5,sprintf('%.2f',max(plot_correlation)),adj=0)
            # text(-8.5,-10,sprintf('%.2f',0),adj=0)
            par(mar=c(5,4,4,2)+0.1)
            dev.off()
        }
    }else if(data_index==2){
        for(pp in seq(10,90,10)){
            turning_percentile = pp
            COL = colorRampPalette(c(colorRampPalette(c('#DDDDDD','#EEEEEE'))(turning_percentile),
                                     colorRampPalette(c('#EEEEEE','#E377C2'))(100-turning_percentile)))(100)
            
            png(paste0(prefix,'Ptr_CellRanger_UMAP_LCM_tricelltype_',pp,'.png'),
                pointsize=10,width=20,height=15,units='cm',res=300)
            par(mar=c(0,0,0,0))
            plot_correlation = correlation_matrix[,'tri_celltype'] %>% as.numeric
            COL_ind = round(plot_correlation/max(plot_correlation)*99+1)
            order_ind = order(plot_correlation)
            plot(Ptr_projection_umap$UMAP.1[order_ind],
                 Ptr_projection_umap$UMAP.2[order_ind],
                 col=COL[COL_ind][order_ind],
                 pch=20,cex=1.4,
                 main='',xlab='UMAP_1',ylab='UMAP_2')
            for(i in seq_along(COL)){
                x = seq(-9,-2,length.out=length(COL)+1)[c(i,i+1)]
                lines(x,rep(-9,2),lwd=40,lend='butt',col=COL[i])
                # y = seq(-10,-5,length.out=length(COL)+1)[c(i,i+1)]
                # lines(rep(-9,2),y,lwd=20,lend='butt',col=COL[i])
                # if(i==pp1) lines(rep(-9.2,2),y,lwd=10,lend='butt')
                # if(i==pp2) lines(rep(-9.2,2),y,lwd=10,lend='butt')
            }
            text(-2,-7.5,sprintf('%.2f',max(plot_correlation)),cex=4)
            text(-9,-7.5,sprintf('%.2f',0),cex=4)
            # text(-9,-4.6,'Cor')
            # text(-8.5,-5,sprintf('%.2f',max(plot_correlation)),adj=0)
            # text(-8.5,-10,sprintf('%.2f',0),adj=0)
            par(mar=c(5,4,4,2)+0.1)
            dev.off()
        }
    }
}




