#install.packages("devtools")
library(devtools) 
#install_github("kassambara/easyGgplot2")
library(easyGgplot2)
#install.packages("reshape2")
#可以把寬數據變成長數據 
library(reshape2)

setwd('E:\\Citrus\\Academy\\Lab YCLin\\Bioinfo')
Ptr_SCseq_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')
#Ptr_SCseq_norm_UMI_matrix = apply(Ptr_SCseq_UMI_matrix,2,function(x)x/sum(x)*1000)
#Ptr_SCseq_log2_norm_UMI_matrix = log2(Ptr_SCseq_norm_UMI_matrix + 1)
Ptr_SCseq_log2_norm_UMI_matrix = readRDS('RDS_Ptr_SCseq_log2_norm_UMI_matrix.rds')

#dim(Ptr_SCseq_UMI_matrix) # 34699 4705
#dim(Ptr_SCseq_norm_UMI_matrix)  # 34699 4705

ANT.sum = Ptr_SCseq_log2_norm_UMI_matrix[c('Potri.001G169500.v4.1','Potri.002G114800.v4.1','Potri.003G064700.v4.1',
                                  'Potri.005G148400.v4.1','Potri.007G007400.v4.1','Potri.014G012200.v4.1'),]
ANT.sum.t = t(ANT.sum) #轉置
#dim(ANT.sum.t)

AtHB8.sum = Ptr_SCseq_log2_norm_UMI_matrix[c('Potri.001G188800.v4.1','Potri.001G372300.v4.1',
                                             'Potri.003G050100.v4.1','Potri.004G211300.v4.1',
                                             'Potri.006G237500.v4.1','Potri.009G014500.v4.1',
                                             'Potri.011G098300.v4.1','Potri.018G045100.v4.1'),]
AtHB8.sum.t = t(AtHB8.sum)#轉置
#dim(AtHB8.sum.t)

VND.sum = Ptr_SCseq_log2_norm_UMI_matrix[c('Potri.015G127400.v4.1','Potri.012G126500.v4.1','Potri.003G113000.v4.1',
                                           'Potri.001G120000.v4.1','Potri.007G014400.v4.1','Potri.005G116800.v4.1'),]
VND.sum.t = t(VND.sum) #轉置

Ray.sum = Ptr_SCseq_log2_norm_UMI_matrix[c('Potri.005G239300.v4.1',
                                           'Potri.013G148900.v4.1',
                                           'Potri.014G162800.v4.1',
                                           'Potri.014G172400.v4.1',
                                           'Potri.002G016000.v4.1',
                                           'Potri.018G052200.v4.1',
                                           'Potri.007G033700.v4.1',
                                           'Potri.014G019500.v4.1',
                                           'Potri.003G111400.v4.1',
                                           'Potri.001G416500.v4.1',
                                           'Potri.002G223100.v4.1',
                                           'Potri.014G106800.v4.1',
                                           'Potri.006G060200.v4.1',
                                           'Potri.008G031700.v4.1',
                                           'Potri.017G114600.v4.1',
                                           'Potri.001G345300.v4.1',
                                           'Potri.001G066400.v4.1',
                                           'Potri.010G247500.v4.1',
                                           'Potri.002G119400.v4.1',
                                           'Potri.013G142760.v4.1',
                                           'Potri.004G162600.v4.1',
                                           'Potri.002G049600.v4.1',
                                           'Potri.005G195600.v4.1',
                                           'Potri.008G065600.v4.1',
                                           'Potri.013G066600.v4.1'
                                           ),]
Ray.sum.t = t(Ray.sum) 

EXPA1.sum = Ptr_SCseq_log2_norm_UMI_matrix[c('Potri.001G240900.v4.1'),]
EXPA1.sum.t = t(EXPA1.sum)
head(EXPA1.sum)

PEAR1.sum = Ptr_SCseq_log2_norm_UMI_matrix[c('Potri.006G084200.v4.1'),]


write.csv(PEAR1.sum,file="PEAR1_log2_norm_UMI.csv",row.names = T)

setwd('E:\\Citrus\\Academy\\Lab YCLin\\SingleCellSeq\\04Analysis\\Expression\\Marker genes')
input_filename="PEAR1_log2_norm_UMI.csv"
data<-read.csv(input_filename, header = T)
df<-melt(data,id.vars = c("Cluster"))
head(df)

###OUTPUT = 900 x 220 px
p <- ggplot2.violinplot(data=df, xName='Cluster',yName='value',
                   xShowTitle=FALSE, ytitle="log2 normalized UMI",
                   xtickLabelRotation=45,xShowTickLabel=FALSE,
                   backgroundColor="white",
                   removePanelGrid=TRUE,removePanelBorder=TRUE,
                   axisLine=c(0.5, "solid", "black"),
                   groupName='Cluster', 
                   groupColors=c("#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A",
                                 "#F8E71C", "#9467BD", "#D62728", "#E377C2",
                                 "#9B9B9B", "#4B4B4B"),
                   faceting=T, facetingVarNames="variable")
                   
p = p + geom_violin(scale='width') + theme(legend.key.size = unit(0.5, 'cm'))
p = p + scale_x_discrete(limits=c("#01","#02","#03","#04","#05","#06","#07","#08"))
p = p + scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
p