library(magrittr)
library(DESeq2)


# Define functions =============================================================
build_sample_table <- function(DEAplan) {
    get_group_label <- function(group, label) {
        Group <- paste(
            rep(
                paste0(label, seq_along(group)),
                times = sapply(Group2Sample[group], length)
            ),
            unlist(sapply(Group2Sample[group], seq_along, simplify = FALSE)),
            sep = "_"
        )
        return(Group)
    }
    ALL_GROUP <- c(DEAplan$EXP_GROUP, DEAplan$CTRL_GROUP)
    out <- data.frame(
        Condition = rep(
            ALL_GROUP,
            times = sapply(Group2Sample[ALL_GROUP], length)
        ),
        Group = c(
            get_group_label(DEAplan$EXP_GROUP, "Target"),
            get_group_label(DEAplan$CTRL_GROUP, "NonTarget")
        ),
        Sample = unlist(Group2Sample[ALL_GROUP], use.names = FALSE)
    )
    return(out)
}


# Set parameters ===============================================================
input_gene_count_matrix_csvs <- snakemake@input$gene_count_matrix_csvs
output_DEA_folder <- snakemake@output$DEA_folder
DEAplan <- snakemake@params$DEAplan
Group2Sample <- snakemake@params$Group2Sample


# Implementation ===============================================================
# Create output directory
if (!dir.exists(output_DEA_folder)) {
    dir.create(output_DEA_folder, recursive = TRUE)
}


sampleTable <- build_sample_table(DEAplan)

read_order <- sapply(
    sampleTable$Sample,
    function(x)
        which(grepl(
            paste0(x, "_gene_count_matrix.csv"), input_gene_count_matrix_csvs
        ))
)
input_csvs <- input_gene_count_matrix_csvs[read_order]
readCountsData_1st <- read.csv(input_csvs[1])
ref_gene_id <- readCountsData_1st$gene_id
readCountsData_res <-
    sapply(
        simplify = FALSE, input_csvs[-1],
        function(x) read.csv(x, row.names = "gene_id")[ref_gene_id, ]
    )
readCountsData <-
    do.call(cbind, c(list(readCountsData_1st), readCountsData_res))
colnames(readCountsData) <-
    sub("_gene_count_matrix.csv", "", basename(colnames(readCountsData)))

isPaired <- (DEAplan$IS_PAIRED == "true")

sigSign <- DEAplan$SIGN

LFC <- log2(DEAplan$FC)
FDR <- as.numeric(DEAplan$FDR)


### Preprocess the raw data and setting
temp <- strsplit(sampleTable$Group, split = "_")
Condition <- sapply(temp, extract, 1)
Replication <- sapply(temp, extract, 2)

geneID <- readCountsData$gene_id
countsMatrix <- readCountsData[, -1] %>%
  round(0) %>%
  as.matrix()
rownames(countsMatrix) <- geneID
countsMatrix %<>% extract(rowSums(.) >= length(Condition), sampleTable$Sample)
filteredGeneID <- rownames(countsMatrix)


### Construct the design matrix
if (isPaired) {
  coldata <- data.frame(
    row.names = sampleTable$Sample,
    Replication = factor(Replication),
    Condition = factor(Condition)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = countsMatrix,
    colData = coldata,
    design = ~ Replication + Condition
  )
} else {
  coldata <- data.frame(
    row.names = sampleTable$Sample,
    Condition = factor(Condition)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = countsMatrix,
    colData = coldata,
    design = ~Condition
  )
}


### Implement the differential expression analysis (DEA)
dds <- DESeq(dds)


### Prepare the output data
wholeOutTable <- as.data.frame(counts(dds, normalized = TRUE)[filteredGeneID, ])
# will be concatenated in following steps
colnames(wholeOutTable) <- paste(sampleTable$Condition, Replication, sep = "_")


### Calculate the group mean counts of each gene
for (C in Condition) {
  Out <- wholeOutTable[, grepl(paste0("^", C), Condition)] %>%
    apply(1, mean) %>%
    as.numeric()
  outName <- sampleTable$Condition[grepl(paste0("^", C), sampleTable$Group)][1]
  code <- paste0("wholeOutTable$mean_", outName, " = Out")
  eval(parse(text = code))
}


### Extract the differentially expressed genes (DEGs)
TargetCondition <- unique(grep("^Target", Condition, value = TRUE))
NonTargetCondition <- setdiff(unique(Condition), TargetCondition)
for (TC in TargetCondition) {
  specificUpGeneID <- filteredGeneID
  specificDownGeneID <- filteredGeneID
  for (NTC in NonTargetCondition) {
    Out <- results(dds,
      alpha = FDR,
      lfcThreshold = 0,
      contrast = c("Condition", TC, NTC)
    )[filteredGeneID, ]

    outTCName <-
      sampleTable$Condition[grepl(paste0("^", TC), sampleTable$Group)][1]
    outNTCName <-
      sampleTable$Condition[grepl(paste0("^", NTC), sampleTable$Group)][1]

    code <- paste0(
      "wholeOutTable$log2FoldChange_", outTCName, "_", outNTCName,
      " = Out$log2FoldChange"
    )
    eval(parse(text = code))

    code <- paste0(
      "wholeOutTable$padj_", outTCName, "_", outNTCName, " = Out$padj"
    )
    eval(parse(text = code))

    temp <- rownames(Out)[with(Out, padj < FDR & log2FoldChange >= LFC)]
    specificUpGeneID %<>% intersect(temp)

    temp <- rownames(Out)[with(Out, padj < FDR & log2FoldChange <= (-LFC))]
    specificDownGeneID %<>% intersect(temp)
  }

  if (sigSign == "both") {
    specificGeneID <- c(specificUpGeneID, specificDownGeneID)
  } else if (sigSign == "up") {
    specificGeneID <- specificUpGeneID
  } else if (sigSign == "down") {
    specificGeneID <- specificDownGeneID
  }

  code <- paste0("specificGeneID_", TC, " = specificGeneID")
  eval(parse(text = code))
}


### Output the whole DEA results table
write.csv(
    wholeOutTable,
    file = file.path(output_DEA_folder, "wholeTable.csv")
)


### Output the DEGs table of each target group
for (TC in TargetCondition) {
  code <- paste0("specificGeneID = specificGeneID_", TC)
  eval(parse(text = code))

  outName <- sampleTable$Condition[grepl(paste0("^", TC), sampleTable$Group)][1]

  input_table <- wholeOutTable[specificGeneID, ]
  isFCcolumn <- grepl(paste0("log2FoldChange_", outName), colnames(input_table))
  reorder_ID <- input_table[, isFCcolumn, drop = FALSE] %>%
    apply(1, sum) %>%
    order(decreasing = TRUE)
  Out <- input_table[reorder_ID, ]

  write.csv(
    Out,
    file = file.path(output_DEA_folder, paste0("DEGtable_", outName, ".csv"))
  )
}
