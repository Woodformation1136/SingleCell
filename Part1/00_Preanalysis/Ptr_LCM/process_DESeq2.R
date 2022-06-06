library(magrittr)
library(DESeq2)

### Input the raw data and setting
Args = commandArgs(trailingOnly=T)

sampleTable = read.csv(Args[1])

readCountsData = read.csv(Args[2])
colnames(readCountsData) %<>% sub('outputStringtie_gtf_','',.)

isPaired = as.logical(Args[3])

sigSign = Args[4]

FDR = as.numeric(Args[5])
LFC = log2(as.numeric(Args[6])) # log2(1.5)


### Preprocess the raw data and setting
temp = strsplit(sampleTable$Group,split='_')
Condition = sapply(temp,extract,1)
Replication =  sapply(temp,extract,2)

geneID = readCountsData$gene_id
countsMatrix = readCountsData[,-1] %>% round(0) %>% as.matrix
rownames(countsMatrix) = geneID
countsMatrix %<>% extract(rowSums(.)>=length(Condition),sampleTable$Group)
filteredGeneID = rownames(countsMatrix)


### Construct the design matrix
if(isPaired){
  coldata = data.frame(row.names=sampleTable$ID,
                       Replication=factor(Replication),
                       Condition=factor(Condition))
  dds = DESeqDataSetFromMatrix(countData=countsMatrix,
                               colData=coldata,
                               design=~Replication+Condition)
}else{
  coldata = data.frame(row.names=sampleTable$ID,
                       Condition=factor(Condition))
  dds = DESeqDataSetFromMatrix(countData=countsMatrix,
                               colData=coldata,
                               design=~Condition)
}


### Implement the differential expression analysis (DEA)
dds = DESeq(dds)


### Prepare the output data
wholeOutTable = as.data.frame(counts(dds,normalized=T)[filteredGeneID,]) # will be concatenate in following steps
colnames(wholeOutTable) = paste(sampleTable$Condition,Replication,sep='_')


### Calculate the group mean counts of each gene
for(C in Condition){
  Out = wholeOutTable[,grepl(paste0('^',C),Condition)] %>% apply(1,mean) %>% as.numeric
  outName = sampleTable$Condition[grepl(paste0('^',C),sampleTable$Group)][1]
  code = paste0("wholeOutTable$mean_",outName," = Out")
  eval(parse(text=code))
}


### Extract the differentially expressed genes (DEGs)
TargetCondition = unique(grep('^Target',Condition,value=T))
NonTargetCondition = setdiff(unique(Condition),TargetCondition)
for(TC in TargetCondition){
  specificUpGeneID = filteredGeneID
  specificDownGeneID = filteredGeneID
  for(NTC in NonTargetCondition){
    Out = results(dds,
                  alpha=FDR,
                  lfcThreshold=0,
                  contrast=c('Condition',TC,NTC))[filteredGeneID,]
    # contrast = c(factor,numerator,denominator)
    
    outTCName = sampleTable$Condition[grepl(paste0('^',TC),sampleTable$Group)][1]
    outNTCName = sampleTable$Condition[grepl(paste0('^',NTC),sampleTable$Group)][1]

    code = paste0("wholeOutTable$log2FoldChange_",outTCName,"_",outNTCName," = Out$log2FoldChange")
    eval(parse(text=code))
    
    code = paste0("wholeOutTable$padj_",outTCName,"_",outNTCName," = Out$padj")
    eval(parse(text=code))
    
    temp = rownames(Out)[with(Out, padj<FDR & log2FoldChange>=LFC)]
    specificUpGeneID %<>% intersect(temp)
    
    temp = rownames(Out)[with(Out, padj<FDR & log2FoldChange<=(-LFC))]
    specificDownGeneID %<>% intersect(temp)
    
  }
  
  if(sigSign=='Both'){
    specificGeneID = c(specificUpGeneID,specificDownGeneID)
  }else if(sigSign=='Up'){
    specificGeneID = specificUpGeneID
  }else if(sigSign=='Down'){
    specificGeneID = specificDownGeneID
  }
  
  code = paste0("specificGeneID_",TC," = specificGeneID")
  eval(parse(text=code))
}


### Output the whole DEA results table
write.csv(wholeOutTable,
          file='outputDESeq2_wholeTable.csv')


### Output the DEGs table of each target group
for(TC in TargetCondition){
  code = paste0("specificGeneID = specificGeneID_",TC)
  eval(parse(text=code))
  
  outName = sampleTable$Condition[grepl(paste0('^',TC),sampleTable$Group)][1]
  
  input_table = wholeOutTable[specificGeneID,]
  isFCcolumn = grepl(paste0('log2FoldChange_',outName),colnames(input_table))
  reorder_ID = input_table[,isFCcolumn,drop=F] %>% 
                apply(1,sum) %>% order(decreasing=T)
  Out = input_table[reorder_ID,]
  
  write.csv(Out,file=paste0('outputDESeq2_DEGtable_',outName,'.csv'))
}

