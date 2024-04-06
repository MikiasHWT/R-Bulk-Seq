## Geo Dataset download and analysis draft script
# setwd("C:/Users/Owner/Desktop/CODE/DB.BRI.BioInfo/March_GeoData/GeoData")

# Load libraries
library(GEOquery)
library(tidyverse)
library(janitor)
library(DESeq2)

# View GEOquery reference manual
# browseVignettes("GEOquery")
# ??GEOquerry


# Load Gene counts

geneCounts <- read.csv("Data/GSE183973_bulkRNA_gene_counts.csv", 
                       check.names = FALSE)

# colnames(geneCounts)[1] <- "genes"

# cleanCounts <- read.csv("Data/GSE183973_bulkRNA_gene_counts.csv", row.names = 1)


################# Load Metadata from Geo
# gse <- getGEO(GEO = "GSE183973", GSEMatrix = TRUE)
# 
# gse
# 
# metadata <- pData(phenoData(gse[[1]]))
# 
# view(metadata)
# colnames(metadata)
# 
# subMetadata <- select(metadata, c("title", "characteristics_ch1", "patient group:ch1"))


# Load Metadata from CSV

 metadata <- read.csv("Data/GSE183973_metadata_samples.csv", check.names = FALSE)



# Clean metadata
# subMetadata <- subMetadata |> 
#   clean_names() |> 
#   mutate(cell_type = str_remove(characteristics_ch1, "cell type: Sorted "), 
#          patient_type = patient_group_ch1) |> 
#   separate(cell_type, into = c("cell_type", "tissue_type"), sep = " from ") |> 
#   select(-characteristics_ch1, -patient_group_ch1)


###################################### TEMPORARY
# # Simplify geneCounts 
# 
# geneCounts <- geneCounts |> 
#   filter(rowSums(across(-X)) > 100)

######################################



# Pivot geneCounts to Long format & Remove leading "X" value in sample strings
 
colnames(geneCounts)[1] <- "genes"

colnames(metadata)[1] <- "samples"


row.names(geneCounts) <- geneCounts$genes

Counts <- geneCounts |>
  select(-genes)


countsLong <- geneCounts |> 
  pivot_longer(cols = !genes, names_to = "samples", values_to = "counts")


featureMatrix <- countsLong |> 
  left_join(metadata, by = c("samples" = "samples"))


# Join Metadata and Counts & reorganize

# Get the desired column order based on metadata$samples
desired_order <- metadata$samples

# Reorder the columns of Counts based on desired_order
Counts <- Counts[, desired_order]




# Finding mismatched names v columns
all(colnames(Counts) %in% metadata$samples)


# Correct order?
all(colnames(Counts) == metadata$samples)



### Exploratory data analysis

# barplot
featureMatrix |> 
  filter(genes == "CD101") |> 
  ggplot(aes(x = samples, y = counts, fill = cell_type)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# density plot
featureMatrix |> 
  filter(genes == "FTL") |> 
  ggplot(aes(x = counts, fill = cell_type)) +
  geom_density(alpha = 0.5)

# boxplot
featureMatrix |> 
  filter(genes == "LYZ") |> 
  ggplot(aes(x = cell_type, y = counts)) +
  geom_boxplot()

# voilinplot 
featureMatrix |> 
  filter(genes == "CD74") |> 
  ggplot(aes(x = cell_type, y = counts)) +
  geom_violin()

# scatterplot
featureMatrix |> 
  filter(genes == "FTL" | genes == "LYZ") |> 
  pivot_wider(names_from = genes, values_from = counts) |> 
  ggplot(aes(x = FTL, y = LYZ, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# heatmap
genesOfInterest <- c("FTL", "FN1", "CD74", "LYZ")

featureMatrix |>  
  filter(genes %in% genesOfInterest) |> 
  ggplot(aes(x = samples, y = genes, fill = counts)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Save plots
# ggsave(plotName, filename = "fileName.pdf", width = 10, height = 8)

## alternatively
# pdf("fileName.pdf, width = 10, height = 8)
# ggplot(BLah Blah Blah)
# dev.off()



######################################################################################
#################FUCK ME SIDEWAYS, Data organization !!!!!!!!! DEFINITLY REDO THIS

# sampleNames <- as.data.frame(colnames(cleanCounts))
# 
# sampleNames <- sampleNames |> 
#   mutate(colNames = str_remove(`colnames(cleanCounts)`, "X")) |> 
#   select(-`colnames(cleanCounts)`)
# 
# colnames(cleanCounts) <- sampleNames$colNames
# 
# arrangedCounts <- cleanCounts |> 
#   select(sort(colnames(cleanCounts)))
#   
# arrangedMetadata <- subMetadata |> 
#   arrange(subMetadata$title)
# 
# rownames(arrangedMetadata) <- arrangedMetadata$title
# 
# arrangedMetadata <- arrangedMetadata  |> 
#   mutate(across(everything(), ~str_replace_all(., "\\s+", ""))) |> 
#   mutate(across(everything(), ~str_replace_all(., "\\+", "pos")))
# 
# 
# # Finding mismatched names v columns
# all(colnames(arrangedCounts) %in% row.names(arrangedMetadata))
# 
# # Correct order?
# all(colnames(arrangedCounts) == row.names(arrangedMetadata))
# 
# arrangedMetadata$patient_type <- factor(arrangedMetadata$patient_type)

###################################################################################### 


dds <- DESeqDataSetFromMatrix(countData = Counts, 
                       colData = metadata, 
                       design = ~ cell_type + patient_group)

dds


# Prefiltering : remove low expressed genes

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]


# set factor level 
## Note: COLLAPSE TECHNICAL REPLICATES
dds$patient_group <- relevel(dds$patient_group, ref = "non_smoker")


# set design criteria & Interactions
# design(dds) <- ~patient_type + tissue_type + patient_type:tissue_type


# run DESeq
dds <- DESeq(dds)

res <- results(dds) # , contrast = c("patient_group", "Non_smoker", "Smoker"))

resultsNames(dds)

head(coef(dds))
res

plotMA(res)

summary(res)


# resSmall <- results(dds, alpha = 0.01)
# 
# summary(resSmall)

# Extract differentially expressed genes
DE_genes <- subset(res, padj < 0.01)

summary(DE_genes)

# MA plot
plotMA(DE_genes)

DEG <- DE_genes@rownames

featureMatrix |>  
  filter(genes %in% DEG) |> 
  ggplot(aes(x = samples, y = genes, fill = counts)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Transform Values for EDA

vsd <- vst(dds, blind=FALSE)

rld <- rlog(dds, blind=FALSE)

ntd <- normTransform(dds)

library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

library("pheatmap")

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)[,c("patient_group","cell_type")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(vsd$cell_type, vsd$patient_group, sep="-")

colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



plotPCA(vsd, intgroup=c("cell_type", "patient_group"))

plotPCA(vsd, intgroup=c("cell_type"))


pcaData <- plotPCA(vsd, 
                   intgroup=c("cell_type", "patient_group"), 
                   returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=patient_group, shape=cell_type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



