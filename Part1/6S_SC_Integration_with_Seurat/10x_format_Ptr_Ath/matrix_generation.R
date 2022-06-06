library(magrittr)
library(dplyr)
library(Matrix)

setwd('C:\\Users\\Jung-Chen\\Desktop\\20210526_10x_format_Ptr_Ath')


#Prepare the ortholog cluster table
ortho_cluster_table = readRDS('RDS_ortho_cluster_table.rds')

Ptr_ortho_cluster_table = filter(ortho_cluster_table,speciesName=='PoT')
Ath_ortho_cluster_table = filter(ortho_cluster_table,speciesName=='ArT')

##The Ath pepetide sequences used for ortholog clustering contains gene isoforms
Ath_geneID = Ath_ortho_cluster_table$geneName %>%
  strsplit(split='[.]') %>% sapply(extract,1)
length(Ath_geneID) # 48359
Ath_isoformID = Ath_ortho_cluster_table$geneName %>%
  strsplit(split='[.]') %>% sapply(extract,2)
length(Ath_isoformID) # 48359

Ath_unqiue_geneID = Ath_ortho_cluster_table$geneName %>%
  strsplit(split='[.]') %>% sapply(extract,1) %>% unique
length(Ath_unqiue_geneID) # 27655

##List 50 genes which have isoforms belonging to different ortholog clusters
# count = 0
# for(gid in seq_along(Ath_unqiue_geneID)){
#   foo = with(Ath_ortho_cluster_table,
#              clusterID[grepl(Ath_unqiue_geneID[gid],geneName)]) %>% unique
#   if(length(foo)>1){
#     count %<>% add(1)
#     print(paste(count,gid,sep=':'))
#     print(filter(Ath_ortho_cluster_table,
#                  grepl(Ath_unqiue_geneID[gid],geneName)))
#   }
#   if(count==50) break
# }

##Check for those gene do not have ".1" isoform
# sum(Ath_isoformID==1) # 27619
# search_id = match(setdiff(Ath_geneID[Ath_isoformID!='1'],
#                           Ath_geneID[Ath_isoformID=='1']),Ath_unqiue_geneID)
# for(gid in search_id){
#   foo = with(Ath_ortho_cluster_table,
#              clusterID[grepl(Ath_unqiue_geneID[gid],geneName)]) %>% unique
#   if(length(foo)>1){
#     print(gid)
#     print(filter(Ath_ortho_cluster_table,
#                  grepl(Ath_unqiue_geneID[gid],geneName)))
#   }
# }

##We decide to use the smallest number isoform to assign ortholog cluster to each gene
max(as.numeric(Ath_isoformID)) # 27
selected_transcriptID = c()
selected_geneID = c()
for(i in 1:27){
  candidate_geneID = Ath_geneID[Ath_isoformID==as.character(i)]
  
  new_geneID = setdiff(candidate_geneID,selected_geneID)
  if(length(new_geneID)>0){
    selected_geneID %<>% c(new_geneID)
    selected_transcriptID %<>% c(paste0(new_geneID,'.',as.character(i)))
  }
  
  if(all(Ath_unqiue_geneID %in% selected_geneID)){
    print(i)
    break
  }
}
stopifnot(length(selected_transcriptID) == length(Ath_unqiue_geneID))
stopifnot(all(Ath_unqiue_geneID %in% 
                sapply(strsplit(selected_transcriptID,'[.]'),extract,1)))

Ath_ortho_cluster_table %<>% filter(geneName%in%selected_transcriptID)
Ath_ortho_cluster_table$transcriptName = Ath_ortho_cluster_table$geneName
Ath_ortho_cluster_table$geneName %<>% strsplit('[.]') %>% sapply(extract,1)



#Read the UMI matrix
##Ptr
Ptr_SCseq_UMI_matrix = readRDS('RDS_Ptr_SCseq_UMI_matrix.rds')


##Ath: 2021_Lopez-Anido
PATH = '2021_Lopez-Anido\\Leaf-HVYNNBBXX\\filtered_gene_bc_matrices\\Arabidopsis\\'
Ath_SCseq_UMI_matrix = readMM(paste0(PATH,'matrix.mtx'))
barcode.names = read.delim(paste0(PATH,'barcodes.tsv'),header=F,stringsAsFactors=F)
feature.names = read.delim(paste0(PATH,'genes.tsv'),header=F,stringsAsFactors=F)
UMI_matrix_geneID = feature.names$V1
colnames(Ath_SCseq_UMI_matrix) = barcode.names$V1
rownames(Ath_SCseq_UMI_matrix) = feature.names$V1

length(UMI_matrix_geneID) # 31110
length(Ath_unqiue_geneID) # 27655
length(intersect(UMI_matrix_geneID,Ath_unqiue_geneID)) # 27628
#The geneID of UMI matrix is inconsistent with ortholog table

# foo = Ath_SCseq_UMI_matrix[!(UMI_matrix_geneID %in% Ath_unqiue_geneID),]
# plot(apply(Ath_SCseq_UMI_matrix,1,mean),
#      apply(Ath_SCseq_UMI_matrix,1,sd))
# points(apply(foo,1,mean),
#        apply(foo,1,sd),col='red')
#Those genes disappear in ortholog table tends to have lower mean and sd
#We decide to only use the geneID overlap between UMI matrix and ortholog table
reshape_and_output = function(Ath_SCseq_UMI_matrix,UMI_matrix_geneID){
  addition_geneID = setdiff(Ath_unqiue_geneID,UMI_matrix_geneID)
  addition_matrix = matrix(0,nrow=length(addition_geneID),ncol=ncol(Ath_SCseq_UMI_matrix),
                           dimnames=list(addition_geneID,colnames(Ath_SCseq_UMI_matrix)))
  out_SCseq_UMI_matrix = rbind(Ath_SCseq_UMI_matrix[UMI_matrix_geneID %in% Ath_unqiue_geneID,],
                               addition_matrix)[Ath_unqiue_geneID,]
  return(as.matrix(out_SCseq_UMI_matrix))
}

Ath_2021L_SCseq_UMI_matrix = reshape_and_output(Ath_SCseq_UMI_matrix,UMI_matrix_geneID)
dim(Ath_2021L_SCseq_UMI_matrix) # 27655 5021
saveRDS(Ath_2021L_SCseq_UMI_matrix,
        'RDS_Ath_2021L_SCseq_UMI_matrix.rds')


##Ath: 2019_Ryu
###Repilcation 1
PATH = '2019_Ryu\\Sample_WT-WERGFP\\filtered_gene_bc_matrices\\TAIR10\\'
Ath_SCseq_UMI_matrix_1 = readMM(paste0(PATH,'matrix.mtx'))
barcode.names = read.delim(paste0(PATH,'barcodes.tsv'),header=F,stringsAsFactors=F)
feature.names = read.delim(paste0(PATH,'genes.tsv'),header=F,stringsAsFactors=F)
UMI_matrix_geneID_1 = feature.names$V1
colnames(Ath_SCseq_UMI_matrix_1) = barcode.names$V1
rownames(Ath_SCseq_UMI_matrix_1) = feature.names$V1
dim(Ath_SCseq_UMI_matrix_1) # 27416 4406

Ath_2019R1_SCseq_UMI_matrix = reshape_and_output(Ath_SCseq_UMI_matrix_1,UMI_matrix_geneID_1)
dim(Ath_2019R1_SCseq_UMI_matrix) # 27655 4406
saveRDS(Ath_2019R1_SCseq_UMI_matrix,
        'RDS_Ath_2019R1_SCseq_UMI_matrix.rds')

###Repilcation 2
PATH = '2019_Ryu\\Sample_WT-WERGFP_2\\filtered_gene_bc_matrices\\TAIR10\\'
Ath_SCseq_UMI_matrix_2 = readMM(paste0(PATH,'matrix.mtx'))
barcode.names = read.delim(paste0(PATH,'barcodes.tsv'),header=F,stringsAsFactors=F)
feature.names = read.delim(paste0(PATH,'genes.tsv'),header=F,stringsAsFactors=F)
UMI_matrix_geneID_2 = feature.names$V1
colnames(Ath_SCseq_UMI_matrix_2) = barcode.names$V1
rownames(Ath_SCseq_UMI_matrix_2) = feature.names$V1
dim(Ath_SCseq_UMI_matrix_2) # 27416 1473

Ath_2019R2_SCseq_UMI_matrix = reshape_and_output(Ath_SCseq_UMI_matrix_2,UMI_matrix_geneID_2)
dim(Ath_2019R2_SCseq_UMI_matrix) # 27655 1473
saveRDS(Ath_2019R2_SCseq_UMI_matrix,
        'RDS_Ath_2019R2_SCseq_UMI_matrix.rds')

###Repilcation 3
PATH = '2019_Ryu\\Sample_WT-WERGFP_3\\filtered_gene_bc_matrices\\TAIR10\\'
Ath_SCseq_UMI_matrix_3 = readMM(paste0(PATH,'matrix.mtx'))
barcode.names = read.delim(paste0(PATH,'barcodes.tsv'),header=F,stringsAsFactors=F)
feature.names = read.delim(paste0(PATH,'genes.tsv'),header=F,stringsAsFactors=F)
UMI_matrix_geneID_3 = feature.names$V1
colnames(Ath_SCseq_UMI_matrix_3) = barcode.names$V1
rownames(Ath_SCseq_UMI_matrix_3) = feature.names$V1
dim(Ath_SCseq_UMI_matrix_3) # 27416 1643

Ath_2019R3_SCseq_UMI_matrix = reshape_and_output(Ath_SCseq_UMI_matrix_3,UMI_matrix_geneID_3)
dim(Ath_2019R3_SCseq_UMI_matrix) # 27655 1643
saveRDS(Ath_2019R3_SCseq_UMI_matrix,
        'RDS_Ath_2019R3_SCseq_UMI_matrix.rds')

all(UMI_matrix_geneID_1==UMI_matrix_geneID_2)
all(UMI_matrix_geneID_1==UMI_matrix_geneID_3)



#Calculate the ortholog cluster expression
overlap_cluster = intersect(unique(Ptr_ortho_cluster_table$clusterID),
                            unique(Ath_ortho_cluster_table$clusterID))
length(overlap_cluster) # 11143

calculate_ortho_UMI_matrix = function(SCseq_UMI_matrix = Ptr_SCseq_UMI_matrix,
                                      ortho_cluster_table = Ptr_ortho_cluster_table,
                                      cal_cluster = overlap_cluster,
                                      output_prefix='Ptr_',
                                      save_RDS=F){
  geneID_list_of_each_cluster = sapply(cal_cluster,
                                       function(cc) with(ortho_cluster_table,
                                                         geneName[clusterID==cc]))
  matrix_geneID = rownames(SCseq_UMI_matrix)
  cluster_indicator_matrix = Matrix(0,
                                    nrow=length(cal_cluster),
                                    ncol=length(matrix_geneID),
                                    dimnames=list(paste0('Cluster_',cal_cluster),
                                                  matrix_geneID))
  for(i in seq_along(geneID_list_of_each_cluster)){
    selected_genes = geneID_list_of_each_cluster[[i]]
    cluster_indicator_matrix[i,selected_genes] = 1
  }
  SCseq_ortho_cluster_UMI_matrix = as.matrix(cluster_indicator_matrix%*%Matrix(SCseq_UMI_matrix))
  if(save_RDS) saveRDS(SCseq_ortho_cluster_UMI_matrix,
                       paste0('RDS_',output_prefix,'SCseq_ortho_cluster_UMI_matrix.rds'))
  return(SCseq_ortho_cluster_UMI_matrix)
}


##Ptr
Ptr_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Ptr_SCseq_UMI_matrix,
                             ortho_cluster_table = Ptr_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Ptr_',
                             save_RDS=F)
dim(Ptr_SCseq_ortho_cluster_UMI_matrix) # 11143 4705


##Ath: 2021_Lopez-Anido
Ath_2021L_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Ath_2021L_SCseq_UMI_matrix,
                             ortho_cluster_table = Ath_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Ath_2021L_',
                             save_RDS=F)
dim(Ath_2021L_SCseq_ortho_cluster_UMI_matrix) # 11143 5021


##Ath: 2019_Ryu
Ath_2019R1_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Ath_2019R1_SCseq_UMI_matrix,
                             ortho_cluster_table = Ath_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Ath_2019R1_',
                             save_RDS=F)
dim(Ath_2019R1_SCseq_ortho_cluster_UMI_matrix) # 11143 4406

Ath_2019R2_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Ath_2019R2_SCseq_UMI_matrix,
                             ortho_cluster_table = Ath_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Ath_2019R2_',
                             save_RDS=F)
dim(Ath_2019R2_SCseq_ortho_cluster_UMI_matrix) # 11143 1473

Ath_2019R3_SCseq_ortho_cluster_UMI_matrix =
  calculate_ortho_UMI_matrix(SCseq_UMI_matrix = Ath_2019R3_SCseq_UMI_matrix,
                             ortho_cluster_table = Ath_ortho_cluster_table,
                             cal_cluster = overlap_cluster,
                             output_prefix='Ath_2019R3_',
                             save_RDS=F)
dim(Ath_2019R3_SCseq_ortho_cluster_UMI_matrix) # 11143 1643




#Prepare the 10x format files for Seurat
output_as_10x_format = function(matrix,prefix){
  
  gz_temp = gzfile(paste0(prefix,'features.tsv.gz'),'w')
  foo = rownames(matrix)
  write.table(data.frame(foo,foo,'Gene Expression'),
              gz_temp,
              sep='\t',row.names=F,col.names=F,quote=F)
  close(gz_temp)
  
  gz_temp = gzfile(paste0(prefix,'barcodes.tsv.gz'),'w')
  write.table(data.frame(paste0(prefix,colnames(matrix))),
              gz_temp,
              sep='\t',row.names=F,col.names=F,quote=F)
  close(gz_temp)
  
  writeMM(Matrix(matrix),paste0(prefix,'matrix.mtx'))
  system(paste0('gzip ',paste0(prefix,'matrix.mtx')))
  
  return(0)
}


# Ptr_SCseq_ortho_cluster_UMI_matrix = readRDS('RDS_Ptr_SCseq_ortho_cluster_UMI_matrix.rds')

output_as_10x_format(Ptr_SCseq_ortho_cluster_UMI_matrix,'Ptr_')
output_as_10x_format(Ath_2021L_SCseq_ortho_cluster_UMI_matrix,'Ath_2021L_')
output_as_10x_format(Ath_2019R1_SCseq_ortho_cluster_UMI_matrix,'Ath_2019R1_')
output_as_10x_format(Ath_2019R2_SCseq_ortho_cluster_UMI_matrix,'Ath_2019R2_')
output_as_10x_format(Ath_2019R3_SCseq_ortho_cluster_UMI_matrix,'Ath_2019R3_')


