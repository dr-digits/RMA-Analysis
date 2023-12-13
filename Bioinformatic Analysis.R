#SRA TOOLKIT####################################
cd C:\Users\User\Documents\sratoolkit.3.0.0-win64\bin
fasterq-dump SRR13759076


#KALLISTO#######################################
kallisto index -i transcripts.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto quant -i transcripts.idx -o output --single -l 200 -s 20 SRR13759076.fastq


#MERGING COUNTS DATA############################
result_dir = 'data/' # The directory that includes the extracted folders
count_tr = data.frame()
tpm_tr = data.frame()
for(i in list.dirs(result_dir,recursive=F)){
  temp = read.csv(paste0(i,'/abundance.tsv'),sep='\t',stringsAsFactors = F)
  temp_count = data.frame(temp$est_counts)
  temp_tpm = data.frame(temp$tpm)
  colnames(temp_count) = gsub(paste0(result_dir,"/"), "", i)
  colnames(temp_tpm) = gsub(paste0(result_dir,"/"), "", i)
  if(ncol(count_tr) == 0){
    count_tr = temp_count
    rownames(count_tr) = temp$target_id
    tpm_tr = temp_tpm
    rownames(tpm_tr) = temp$target_id
  } else {
    count_tr = cbind(count_tr, temp_count)
    tpm_tr = cbind(tpm_tr,temp_tpm)
  }
}
biomart_file = 'mart_export.txt' #adjust this based on your file location
mapping = read.csv(biomart_file,sep='\t',stringsAsFactors = F,row.names = 1)
count_gn = merge(count_tr,mapping['Gene.stable.ID'], by=0,all=F) # merge only for the shared row names
count_gn = count_gn[,2:ncol(count_gn)]
count_gn = aggregate(.~Gene.stable.ID, count_gn, sum)
rownames(count_gn) = count_gn$Gene.stable.ID 
tpm_gn = merge(tpm_tr,mapping['Gene.stable.ID'], by=0,all=F) # merge only for the shared row names
tpm_gn = tpm_gn[,2:ncol(tpm_gn)]
tpm_gn = aggregate(.~Gene.stable.ID, tpm_gn, sum)
rownames(tpm_gn) = tpm_gn$Gene.stable.ID
write.table(count_gn,file='/path_to_save/count_gene.txt',sep = '\t', na = '',row.names = F)
write.table(tpm_gn,file='/path_to_save/tpm_gene.txt',sep = '\t', na = '',row.names = F)


#DIFFERENTIAL EXPRESSION ANALYSIS###############
library("tidyverse")
library("DESeq2")
metadata = read.csv("metadata_file.csv",row.names = 1)
count = read.csv("count_file.csv",row.names = 1)
# making sure the row names in colData matches to column names in counts_data
all(colnames(count) %in% rownames(metadata))
# are they in the same order?
all(colnames(count) == rownames(metadata))
conds=as.factor(metadata$Group)
coldata <- data.frame(row.names=rownames(metadata),conds)
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(count)),colData=coldata,design=~conds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
cond1 = 'Sarcopenia' #First Condition
cond2 = 'Healthy' #Reference Condition
res=results(dds,contrast=c('conds',cond1,cond2))
res=data.frame(res)
write.table(res,file='sarcopenia_vs_healthy.txt',sep = '\t', na = '',row.names = T,col.names=NA)


#GENE ONTOLOGY ANALYSIS#########################
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
genelist_o_s <- read_csv("genelist_h_s.csv")
genes <- genelist_h_s$gene
OrgDb <- org.Hs.eg.db
go_enrich <- enrichGO(genes, OrgDb, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05)
write.csv(go_enrich,"results_goa_o_s.csv")


#CO-EXPRESSION NETWORK ANALYSIS#################
library(WGCNA)
library(DESeq2)
library(CorLevelPlot)
library(gridExtra)
library(tidyverse)

data <- read.csv('counts.csv', row.names = 1)
phenoData <- read.csv('metadata.csv', row.names = 1)

all(rownames(phenoData) %in% colnames(data))
all(rownames(phenoData) == colnames(data))

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

data <- data[gsg$goodGenes == TRUE,]

htree <- hclust(dist(t(data)), method = "average")
plot(htree)

samples.to.be.excluded <- c('X_41', 'X_78')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)

all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1)

dds75 <- dds[rowSums(counts(dds) >= 15) >= 38,]
nrow(dds75)

dds_norm <- vst(dds75)

norm.counts <- assay(dds_norm) %>% 
  t()

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 3
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 30000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3,
                          minCoreKME = 0.5)

cor <- temp_cor

module_eigengenes <- bwnet$MEs

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('Sarcopenia', GROUP), 1, 0)) %>% 
  select(2)

nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[20],
             y = names(heatmap.data)[1:19],
             col = c("blue1", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
module.membership.measure.pvals[1:10,1:10]

gene.signf.corr <- cor(norm.counts, traits$disease_state_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

write.csv(table(bwnet$colors), 'modules.csv')
write.csv(module.trait.corr.pvals, 'module.trait.corr.pvals.csv')
write.csv(module.trait.corr, 'module.trait.corr.csv')
write.csv(module.gene.mapping, 'module.gene.mapping.csv')
write.csv(t(module.membership.measure), 'module.membership.measure.csv')


#INTERSECT####################################
library('tidyverse')

upregulated_DEGs <- read_csv("up_degs.csv")
top1_module <- read_csv("1%_module.csv")


commons <- Reduce(intersect, list(upregulated_DEGs$...1, top1_module$node))

x <- tibble(commons)
 
write.csv(x, "commons.csv")