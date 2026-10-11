## Set your working directory

# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)

library(tidyverse)

library(CorLevelPlot) 
library(gridExtra)

library(Rmisc) 


# 1. Import the counts matrix and metadata and filter using DESeq2

countsTable <- read.table("mydata/salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "mydata/Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 



# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# All of our samples look good; the genes that don't look good shouldn't cluster well and 
# should end up in the "grey", catch all module for uncorrelated genes.

# 3. Normalization ----------------------------------------------------------------------

colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order



# perform variance stabilization
dds_norm <- vst(dds)
# dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices 

# visualization to pick power

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
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power_6 <- 6 # should try running with 7 based on R^2
soft_power_4 <- 4 # should try running with 7 based on R^2
soft_power_8 <- 8 # should try running with 7 based on R^2
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package


# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.
#bwnet <- blockwiseModules(norm.counts,
#                          maxBlockSize = 26000,
 #                         minModuleSize = 30, 
  #                        reassignThreshold=0,
   #                       TOMType = "signed",
    #                      power = soft_power,
     #                     mergeCutHeight = 0.25,
      #                    numericLabels = F,
       #               verbose = 3)

# TOMtype (Topological Overlap Matrix type) parameter - unsigned - doesn't consider positive/negative co-expression
# signed - when you want to consider the direction of co-expression interaction, e.g., activating or inhibiting
# WGCNA often uses a dendrogram-based approach to identify modules. The choice of the 
# height cut in the dendrogram can determine the number of modules. Selecting a higher
# cut height results in fewer, larger modules, while a lower cut height leads to more, 
# smaller modules.

saveRDS(bwnet, file = "mydata/bwnet.rds")

bwnet.4 <- readRDS("mydata/bwnet.4.rds")
bwnet.6 <- readRDS("mydata/bwnet.6.rds")
bwnet.8 <- readRDS("mydata/bwnet.8.rds")


cor <- temp_cor
This step above takes a really long time and a lot of memory - this is the creation/assignment of all of the modules. If it doesn’t finish within 3 minutes of class, please copy the file I placed in /gpfs1/cl/biol3990/Transcriptomics/WGCNA called bwnet.rds.



# 5. Module Eigengenes ---------------------------------------------------------


#4
module_eigengenes.4 <- bwnet.4$MEs

head(module_eigengenes.4)


# get number of genes for each module
table(bwnet.4$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet.4$dendrograms[[1]], cbind(bwnet.4$unmergedColors, bwnet.4$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
#6
module_eigengenes.6 <- bwnet.6$MEs

head(module_eigengenes.6)


# get number of genes for each module
table(bwnet.6$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet.6$dendrograms[[1]], cbind(bwnet.6$unmergedColors, bwnet.6$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
#8
module_eigengenes.8 <- bwnet.8$MEs

head(module_eigengenes.8)


# get number of genes for each module
table(bwnet.8$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet.8$dendrograms[[1]], cbind(bwnet.8$unmergedColors, bwnet.8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module


# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

traits <- sample_metadata[, c(5,8,11,14,17)]


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#4
module.trait.corr.4 <- cor(module_eigengenes.4, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr.4, nSamples)
#6
module.trait.corr.6 <- cor(module_eigengenes.6, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr.6, nSamples)
#8
module.trait.corr.8 <- cor(module_eigengenes.8, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr.8, nSamples)



# visualize module-trait association as a heatmap
#4
heatmap.data.4 <- merge(module_eigengenes.4, traits, by = 'row.names')
#6
heatmap.data.6 <- merge(module_eigengenes.6, traits, by = 'row.names')
#8
heatmap.data.8 <- merge(module_eigengenes.8, traits, by = 'row.names')

#4
head(heatmap.data.4)
#6
head(heatmap.data.6)
#8
head(heatmap.data.8)

#4
heatmap.data.4 <- heatmap.data.4 %>% 
  column_to_rownames(var = 'Row.names')
#6
heatmap.data.6 <- heatmap.data.6 %>% 
  column_to_rownames(var = 'Row.names')

#8
heatmap.data.8 <- heatmap.data.8 %>% 
  column_to_rownames(var = 'Row.names')

names(heatmap.data.4)
names(heatmap.data.6)
names(heatmap.data.8)

#4
CorLevelPlot(heatmap.data.4,
             x = names(heatmap.data.4)[11:15],
             y = names(heatmap.data.4)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))
#6
CorLevelPlot(heatmap.data.6,
             x = names(heatmap.data.6)[12:16],
             y = names(heatmap.data.6)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))
#8
CorLevelPlot(heatmap.data.8,
             x = names(heatmap.data.8)[16:20],
             y = names(heatmap.data.8)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))


#4
module.gene.mapping.4 <- as.data.frame(bwnet.4$colors) # assigns module membership to each gene
#6
module.gene.mapping.6 <- as.data.frame(bwnet.6$colors)
#8
module.gene.mapping.8 <- as.data.frame(bwnet.8$colors)


module.gene.mapping.4 %>% 
  filter(`bwnet.4$colors` == 'yellow') %>% 
  rownames()

groups <- sample_metadata[,c(3,1)]
module_eigengene.metadata.4 <- merge(groups, heatmap.data.4, by = 'row.names')

#Create a summary data frame of a particular module eigengene information
MEyellow_summary.4 <- summarySE(module_eigengene.metadata.4, measurevar="MEyellow", groupvars=c("Generation","treatment"))

#Plot a line interaction plot of a particular module eigengene
ggplot(MEyellow_summary.4, aes(x=as.factor(Generation), y=MEyellow, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=MEyellow-se, ymax=MEyellow+se), width=.15) +
  geom_line(aes(color=treatment, group=treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))


# 6B. Intramodular analysis: Identifying driver genes ---------------

# Get top hub genes (genes with highest connectivity in the network)
hubs  <-  chooseTopHubInEachModule(norm.counts, bwnet.4$colors, type = "signed", omitColors = "")
hubs

### Plot Individual genes  to check! ### 

d <-plotCounts(dds, gene="TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434", intgroup = (c("treatment","Generation")), returnData=TRUE)
d_summary <- summarySE(d, measurevar = "count", groupvars=c("Generation","treatment"))

ggplot(d_summary, aes(x=as.factor(Generation), y=count, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=count-se, ymax=count+se), width=.15) +
  geom_line(aes(color=treatment, group=treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))

# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership
library(pheatmap)
t_norm.counts <- norm.counts %>% t() %>% as.data.frame()

# purple module
purple_transcripts <- module.gene.mapping.4 %>% 
  filter(`bwnet.4$colors` == 'purple') %>% 
  rownames()

t_norm.counts_purple <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% purple_transcripts)

t_norm.counts_purple <- t_norm.counts_purple - rowMeans(t_norm.counts_purple)
df <- as.data.frame(colData(dds)[,c("Generation","treatment")])

#blue to yellow heatmap color scheme; plot purple module genes
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaks <- c(seq(min(t_norm.counts_purple), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_purple)/paletteLength, max(t_norm.counts_purple), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_purple, color = myColor, breaks = myBreaks,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "purple")


# yellow module
yellow_transcripts <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()

t_norm.counts_yellow <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% yellow_transcripts)

t_norm.counts_yellow <- t_norm.counts_yellow - rowMeans(t_norm.counts_yellow)
df <- as.data.frame(colData(dds)[,c("generation","treatment")]) # Change generation (values) to Generations (categories)

#blue to yellow heatmap color scheme; yellow module
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaks <- c(seq(min(t_norm.counts_yellow), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_yellow)/paletteLength, max(t_norm.counts_yellow), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_yellow, color = myColor, breaks = myBreaks,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "yellow")





# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure.4 <- cor(module_eigengenes.4, norm.counts, use = 'p')
module.membership.measure.pvals.4 <- corPvalueStudent(module.membership.measure.4, nSamples)


module.membership.measure.pvals.4[1:10,1:10]



###########################################
# Saving for the future and for TopGO analyses
###########################################

module.gene.mapping$geneID <- rownames(module.gene.mapping)

parts <- strsplit(module.gene.mapping$geneID, "::")

module.gene.mapping$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse = "::")
)

write.csv(
  module.gene.mapping,
  "myresults/WGCNA_module_gene_mapping_w_shortIDs.csv",
  row.names = FALSE
)

############################################################
# Extract yellow module genes for TopGO
############################################################

yellow_shortIDs.4 <- unique(
  module.gene.mapping.4$shortID[
    module.gene.mapping.4$`bwnet.4$colors` == "yellow"
  ]
)

length(yellow_shortIDs.4)

head(yellow_shortIDs.4)


write.csv(
  data.frame(shortID = yellow_shortIDs),
  "myresults/yellow_module_genes.csv",
  row.names = FALSE
)