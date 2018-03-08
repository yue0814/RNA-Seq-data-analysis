# 20170516_NA7DTWYE_Gene_Expression
library(SummarizedExperiment)
library(edgeR)
library(DESeq2, quietly = TRUE)
library(biomaRt)
library(systemPipeR)
setwd("I:/a4703/Analysis/4.Clustering/Yue/DESeq")

dat <- readRDS("20170516_NA7DTWYE_Gene_Expression.rds")
col_data <- as.data.frame(colData(dat))
count <- assay(dat, 'exp_count')
countMatrix <- as.matrix(count)



# Split into 2-3-cluster
cluster_data <- read.csv("Cases and clusters from combined DECAMP analysis 20171128.csv")
colnames(cluster_data) <- c("Case", "X2_means_cluster", "X3_means_cluster")
col_data$INDIVIDUAL_ID
#cluster_data <- cluster_data[,c("Case", "X2_means_cluster", "X3_means_cluster")]
cluster_data$Case <- as.character(cluster_data$Case)
study <- c()
case <- c()
for (i in c(1:length(cluster_data$Case))){
  study <- c(study, c(strsplit(cluster_data$Case[i], "-")[[1]])[1])
  case <- c(case, c(strsplit(cluster_data$Case[i], "-")[[1]])[2])
}
cluster_data$Case <- case
cluster_data$Study <- study

for (i in c(1:length(cluster_data$Case))){
  if (length(cluster_data$Case[i])<4){
    cluster_data$Case[i] <- paste0(paste(rep('0',4-nchar(cluster_data$Case[i])),collapse = ''),cluster_data$Case[i])
  }
}
# label the cluster number
col_data$X2_means_cluster <- rep(0,length(col_data$INDIVIDUAL_ID))
col_data$X3_means_cluster <- rep(0,length(col_data$INDIVIDUAL_ID))
for (i in c(1:length(col_data$INDIVIDUAL_ID))){
  cases <- c(strsplit(col_data$INDIVIDUAL_ID[i], "-")[[1]])
  idx <- which(cluster_data$Study==cases[2]&cluster_data$Case==cases[3])
  if (length(idx)>0){
    col_data$X2_means_cluster[i] <- cluster_data$X2_means_cluster[idx]
    col_data$X3_means_cluster[i] <- cluster_data$X3_means_cluster[idx]
  }
}
# remove the row with no cluster membership number
countMatrix <- countMatrix[,-c(which(col_data$X2_means_cluster==0))]
countMatrix <- matrix(as.integer(unlist(countMatrix)), nrow = dim(count)[1])
col_data <- col_data[-c(which(col_data$X2_means_cluster==0)),]
col_data$X2_means_cluster <- as.factor(col_data$X2_means_cluster)
col_data$X3_means_cluster <- as.factor(col_data$X3_means_cluster)


# Sample-wise correlation analysis
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = col_data, design = ~ X2_means_cluster)
d <- cor(assay(rlog(dds)), method="spearman")
save(d, file = "d.RData")
rownames(d) <- c(1:151)
hc <- hclust(dist(1-d))

#pdf("sample_tree.pdf", width=10, height = 30)
#plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
#dev.off()


###### analysis of differentially expressed genes (DEGs) and GO 2017.12.07 ############################
#source("https://bioconductor.org/biocLite.R")
###############DATA MANAGEMENT#######################
targets <- col_data[,c("INDIVIDUAL_ID", "X2_means_cluster")]
row.names(targets) <- NULL
#colnames(targets) <- c("INDIVIDUAL_ID", "Factor")
targets$Factor <- rep(0, 151)
for (i in c(1:length(targets$INDIVIDUAL_ID))){
  targets$Factor[i] <- as.character(targets$X2_means_cluster[i])
}
for (i in c(1:length(targets$INDIVIDUAL_ID))){
  targets$Factor[i] <- paste0("G",as.character(targets$Factor[i]))
}

count_new <- as.data.frame(countMatrix)
rownames(count_new) <- rownames(count)
# push function
push <- function(l, x) {
  lst <- get(l, parent.frame())
  lst[length(lst)+1] <- x
  assign(l, lst, envir=parent.frame())
}
n1 <- list(targets$INDIVIDUAL_ID[1])
for (i in c(2:length(targets$INDIVIDUAL_ID))){
  push('n1', targets$INDIVIDUAL_ID[i])
}
# change "-" to "_"
for (i in c(1:length(targets$INDIVIDUAL_ID))){
  n1[i] <- gsub("-", "_", n1[[i]])
}
names(count_new) <- n1

for (i in c(1:length(targets$INDIVIDUAL_ID))){
  targets$INDIVIDUAL_ID[i] <- gsub("-", "_", targets$INDIVIDUAL_ID[[i]])
}

targets <- targets[, c("INDIVIDUAL_ID", "Factor")]
targets$INDIVIDUAL_ID <- as.factor(targets$INDIVIDUAL_ID)
targets$Factor <- as.factor(targets$Factor)
#write.table(targets, file = "targets.txt", sep = ' ')

# remove duplicated
idx <- which(duplicated(unlist(n1)))
targets <- targets[-c(idx),]
count_new <- count_new[,-c(idx)]
names(targets) <- c("SampleName", "Factor")
tg <- list(CMPSet1 = matrix(c("G2", "G1"), ncol=2))
###############DATA MANAGEMENT#######################

# analysis
#edgeDF <- run_edgeR(countDF=count_new, targets=targets, cmp=tg[[1]],independent=FALSE, mdsplot="")
#DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=10), plot=T)

##############################################################################
#### 2 clusters
degseqDF <- run_DESeq2(countDF=count_new, targets=targets, cmp=tg[[1]], independent=FALSE)

#
# where Fold refers to the fold change cutoff (unlogged) and FDR to the p-value cutoff.
DEG_list2 <- filterDEGs(degDF=degseqDF, filter=c(Fold=1, FDR=10))

# make the table
table1 <- degseqDF[order(degseqDF$`G2-G1_FDR`),]
table1 <- table1[which(table1$`G2-G1_FDR`<0.1),]
table1_greater_0 <- table1[which(table1$`G2-G1_logFC`>0),]
table1_less_0 <- table1[which(table1$`G2-G1_logFC`<=0),]

table1 <- rbind(table1_greater_0, table1_less_0)
#write.csv(table1, file="2_cluster_results.csv")
# design CATDB
ensembl=useMart("ensembl")
listDatasets(ensembl)
# hsapiens_gene_ensembl   Human genes (GRCh38.p10)   
#ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

filters = listFilters(ensembl)
filters[1:5,]

# Attributes define the values we are interested in to retrieve
attributes = listAttributes(ensembl)
attributes[1:5,]

# attributes problem 

#go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
#go <- getBM(attributes=c('go_id', "go_linkage_type" ,'namespace_1003'), mart=ensembl)
go <- getBM(attributes=c("go_id", 'ensembl_gene_id', 'namespace_1003'), mart=ensembl)
goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=ensembl)[,1])

go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
#write.table(go, "GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
catdb <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)


#GO TERM
up_down <- DEG_list2$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")
up <- DEG_list2$Up; names(up) <- paste(names(up), "_up", sep="")
down <- DEG_list2$Down; names(down) <- paste(names(down), "_down", sep="")
DEGlist <- c(up_down, up, down)
DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type="gene", CLSZ=2, cutoff=0.10, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.10, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)


gos <- BatchResultslim[grep("G2-G1_up_down", BatchResultslim$CLID), ]
gos <- BatchResultslim
pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
goBarplot(gos, gocat="BP")
goBarplot(gos, gocat="CC")
dev.off()

##############################################################################
## 3 clusters
targets3 <- col_data[,c("INDIVIDUAL_ID", "X3_means_cluster")]
row.names(targets3) <- NULL
#colnames(targets) <- c("INDIVIDUAL_ID", "Factor")
targets3$Factor <- rep(0, 151)
for (i in c(1:length(targets3$INDIVIDUAL_ID))){
  targets3$Factor[i] <- as.character(targets3$X3_means_cluster[i])
}
for (i in c(1:length(targets3$INDIVIDUAL_ID))){
  targets3$Factor[i] <- paste0("G",as.character(targets3$Factor[i]))
}



for (i in c(1:length(targets3$INDIVIDUAL_ID))){
  targets3$INDIVIDUAL_ID[i] <- gsub("-", "_", targets3$INDIVIDUAL_ID[[i]])
}

targets3 <- targets3[, c("INDIVIDUAL_ID", "Factor")]
targets3$INDIVIDUAL_ID <- as.factor(targets3$INDIVIDUAL_ID)
targets3$Factor <- as.factor(targets3$Factor)

idx <- which(duplicated(unlist(n1)))
targets3 <- targets3[-c(idx),]
names(targets3) <- c("SampleName", "Factor")
tg3 <- list(CMPSet1 = matrix(c("G2", "G3", "G3", "G1", "G1", "G2"), ncol=2))

#edgeDF3 <- run_edgeR(countDF=count_new, targets=targets3, cmp=tg3[[1]],independent=FALSE, mdsplot="")
#DEG3_list <- filterDEGs(degDF=edgeDF3, filter=c(Fold=2, FDR=50), plot=FALSE)
#names(DEG3_list)

degseqDF3 <- run_DESeq2(countDF=count_new, targets=targets3, cmp=tg3[[1]], independent=FALSE)
DEG3_list2 <- filterDEGs(degDF=degseqDF3, filter=c(Fold=1, FDR=10))

table2_1 <- degseqDF3[which(degseqDF3$`G2-G1_FDR`<0.1),c("G2-G1_baseMean", "G2-G1_logFC",
                                                            "G2-G1_lfcSE","G2-G1_stat","G2-G1_pvalue","G2-G1_FDR" )]
table2_1_greater_0 <- table2_1[which(table2_1$`G2-G1_logFC`>0),]
table2_1_greater_0 <- table2_1_greater_0[order(table2_1_greater_0$`G2-G1_FDR`),]
table2_1_less_0 <- table2_1[which(table2_1$`G2-G1_logFC`<=0),]
table2_1_less_0 <- table2_1_less_0[order(table2_1_less_0$`G2-G1_FDR`),]
table2_1 <- rbind(table2_1_greater_0, table2_1_less_0)
#write.csv(table2_1, file="3_cluster_results2_1.csv")
table3_1 <- degseqDF3[which(degseqDF3$`G3-G1_FDR`<0.1),c("G3-G1_baseMean", "G3-G1_logFC",
                                                         "G3-G1_lfcSE","G3-G1_stat","G3-G1_pvalue","G3-G1_FDR" )]
table3_1_greater_0 <- table3_1[which(table3_1$`G3-G1_logFC`>0),]
table3_1_greater_0 <- table3_1_greater_0[order(table3_1_greater_0$`G3-G1_FDR`),]
table3_1_less_0 <- table3_1[which(table3_1$`G3-G1_logFC`<=0),]
table3_1_less_0 <- table3_1_less_0[order(table3_1_less_0$`G3-G1_FDR`),]
table3_1 <- rbind(table3_1_greater_0, table3_1_less_0)
#write.csv(table3_1, file="3_cluster_results3_1.csv")
table3_2 <- degseqDF3[which(degseqDF3$`G3-G2_FDR`<0.1),c("G3-G2_baseMean", "G3-G2_logFC",
                                                         "G3-G2_lfcSE","G3-G2_stat","G3-G2_pvalue","G3-G2_FDR" )]
table3_2_greater_0 <- table3_2[which(table3_2$`G3-G2_logFC`>0),]
table3_2_greater_0 <- table3_2_greater_0[order(table3_2_greater_0$`G3-G2_FDR`),]
table3_2_less_0 <- table3_2[which(table3_2$`G3-G2_logFC`<=0),]
table3_2_less_0 <- table3_2_less_0[order(table3_2_less_0$`G3-G2_FDR`),]
table3_2 <- rbind(table3_2_greater_0, table3_2_less_0)
#write.csv(table3_2, file="3_cluster_results3_2.csv")

up_down <- DEG3_list2$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")
up <- DEG3_list2$Up; names(up) <- paste(names(up), "_up", sep="")
down <- DEG3_list2$Down; names(down) <- paste(names(down), "_down", sep="")
DEGlist <- c(up_down, up, down)
DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type="gene", CLSZ=2, cutoff=0.10, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.10, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)


gos <- BatchResultslim[grep("G3-G2-G1_up_down", BatchResultslim$CLID), ]
gos <- BatchResultslim
pdf("GOslimbarplotMF3.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
goBarplot(gos, gocat="BP")
goBarplot(gos, gocat="CC")
dev.off()
