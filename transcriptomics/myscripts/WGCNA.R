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
#This step above takes a really long time and a lot of memory - this is the creation/assignment of all of the modules. If it doesn’t finish within 3 minutes of class, please copy the file I placed in /gpfs1/cl/biol3990/Transcriptomics/WGCNA called bwnet.rds.



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

####

###4
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

# turquoise module
turquoise_transcripts <- module.gene.mapping.4 %>% 
  filter(`bwnet.4$colors` == 'turquoise') %>% 
  rownames()

t_norm.counts_turquoise <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% turquoise_transcripts)

t_norm.counts_turquoise <- t_norm.counts_turquoise - rowMeans(t_norm.counts_turquoise)
df <- as.data.frame(colData(dds)[,c("Generation","treatment")])

#blue to yellow heatmap color scheme; plot purple module genes
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaksTurquoise <- c(seq(min(t_norm.counts_turquoise), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_turquoise)/paletteLength, max(t_norm.counts_turquoise), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_turquoise, color = myColor, breaks = myBreaksTurquoise,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Turquoise.4")


# yellow module
yellow_transcripts <- module.gene.mapping.4 %>% 
  filter(`bwnet.4$colors` == 'yellow') %>% 
  rownames()

t_norm.counts_yellow <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% yellow_transcripts)

t_norm.counts_yellow <- t_norm.counts_yellow - rowMeans(t_norm.counts_yellow)
df <- as.data.frame(colData(dds)[,c("generation","treatment")]) # Change generation (values) to Generations (categories)

#blue to yellow heatmap color scheme; yellow module
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaksYellow <- c(seq(min(t_norm.counts_yellow), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_yellow)/paletteLength, max(t_norm.counts_yellow), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_yellow, color = myColor, breaks = myBreaksYellow,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Yellow.4")

####

###6
module.gene.mapping.6 %>% 
  filter(`bwnet.6$colors` == 'yellow') %>% 
  rownames()

groups <- sample_metadata[,c(3,1)]
module_eigengene.metadata.6 <- merge(groups, heatmap.data.6, by = 'row.names')

#Create a summary data frame of a particular module eigengene information
MEyellow_summary.6 <- summarySE(module_eigengene.metadata.6, measurevar="MEyellow", groupvars=c("Generation","treatment"))

#Plot a line interaction plot of a particular module eigengene
ggplot(MEyellow_summary.6, aes(x=as.factor(Generation), y=MEyellow, color=treatment, fill = treatment, shape = treatment)) +
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
hubs.6  <-  chooseTopHubInEachModule(norm.counts, bwnet.8$colors, type = "signed", omitColors = "")
hubs.6

### Plot Individual genes  to check! ### 

d.6 <-plotCounts(dds, gene="TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434", intgroup = (c("treatment","Generation")), returnData=TRUE)
d_summary.6 <- summarySE(d.6, measurevar = "count", groupvars=c("Generation","treatment"))

ggplot(d_summary.6, aes(x=as.factor(Generation), y=count, color=treatment, fill = treatment, shape = treatment)) +
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

# turquoise module
turquoise_transcripts.6 <- module.gene.mapping.6 %>% 
  filter(`bwnet.6$colors` == 'turquoise') %>% 
  rownames()

t_norm.counts_turquoise.6 <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% turquoise_transcripts.6)

t_norm.counts_turquoise.6 <- t_norm.counts_turquoise.6 - rowMeans(t_norm.counts_turquoise.6)
df <- as.data.frame(colData(dds)[,c("Generation","treatment")])

#blue to yellow heatmap color scheme; plot purple module genes
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaksTurquoise.6 <- c(seq(min(t_norm.counts_turquoise.6), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(max(t_norm.counts_turquoise.6)/paletteLength, max(t_norm.counts_turquoise.6), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_turquoise.6, color = myColor, breaks = myBreaksTurquoise.6,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Turquoise.6")


# yellow module
yellow_transcripts.6 <- module.gene.mapping.6 %>% 
  filter(`bwnet.6$colors` == 'yellow') %>% 
  rownames()

t_norm.counts_yellow.6 <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% yellow_transcripts.6)

t_norm.counts_yellow.6 <- t_norm.counts_yellow.6 - rowMeans(t_norm.counts_yellow.6)
df <- as.data.frame(colData(dds)[,c("generation","treatment")]) # Change generation (values) to Generations (categories)

#blue to yellow heatmap color scheme; yellow module
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaksYellow.6 <- c(seq(min(t_norm.counts_yellow.6), 0, length.out=ceiling(paletteLength/2) + 1), 
                    seq(max(t_norm.counts_yellow.6)/paletteLength, max(t_norm.counts_yellow.6), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_yellow.6, color = myColor, breaks = myBreaksYellow.6,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Yellow.6")

####

###8
module.gene.mapping.8 %>% 
  filter(`bwnet.8$colors` == 'yellow') %>% 
  rownames()

groups <- sample_metadata[,c(3,1)]
module_eigengene.metadata.8 <- merge(groups, heatmap.data.8, by = 'row.names')

#Create a summary data frame of a particular module eigengene information
MEyellow_summary.8 <- summarySE(module_eigengene.metadata.8, measurevar="MEyellow", groupvars=c("Generation","treatment"))

#Plot a line interaction plot of a particular module eigengene
ggplot(MEyellow_summary.8, aes(x=as.factor(Generation), y=MEyellow, color=treatment, fill = treatment, shape = treatment)) +
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
hubs.8  <-  chooseTopHubInEachModule(norm.counts, bwnet.8$colors, type = "signed", omitColors = "")
hubs.8

### Plot Individual genes  to check! ### 

d.8 <-plotCounts(dds, gene="TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434", intgroup = (c("treatment","Generation")), returnData=TRUE)
d_summary.8 <- summarySE(d.8, measurevar = "count", groupvars=c("Generation","treatment"))

ggplot(d_summary.8, aes(x=as.factor(Generation), y=count, color=treatment, fill = treatment, shape = treatment)) +
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

# turquoise module
turquoise_transcripts.8 <- module.gene.mapping.8 %>% 
  filter(`bwnet.8$colors` == 'turquoise') %>% 
  rownames()

t_norm.counts_turquoise.8 <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% turquoise_transcripts.8)

t_norm.counts_turquoise.8 <- t_norm.counts_turquoise.8 - rowMeans(t_norm.counts_turquoise.8)
df <- as.data.frame(colData(dds)[,c("Generation","treatment")])

#blue to yellow heatmap color scheme; plot purple module genes
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaksTurquoise.8 <- c(seq(min(t_norm.counts_turquoise.8), 0, length.out=ceiling(paletteLength/2) + 1), 
                         seq(max(t_norm.counts_turquoise.8)/paletteLength, max(t_norm.counts_turquoise.8), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_turquoise.8, color = myColor, breaks = myBreaksTurquoise.8,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Turquoise.8")


# yellow module
yellow_transcripts.8 <- module.gene.mapping.8 %>% 
  filter(`bwnet.8$colors` == 'yellow') %>% 
  rownames()

t_norm.counts_yellow.8 <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% yellow_transcripts.8)

t_norm.counts_yellow.8 <- t_norm.counts_yellow.8 - rowMeans(t_norm.counts_yellow.8)
df <- as.data.frame(colData(dds)[,c("generation","treatment")]) # Change generation (values) to Generations (categories)

#blue to yellow heatmap color scheme; yellow module
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaksYellow.8 <- c(seq(min(t_norm.counts_yellow.8), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(t_norm.counts_yellow.8)/paletteLength, max(t_norm.counts_yellow.8), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_yellow.8, color = myColor, breaks = myBreaksYellow.8,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Yellow.8")








# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
#4
module.membership.measure.4 <- cor(module_eigengenes.4, norm.counts, use = 'p')
module.membership.measure.pvals.4 <- corPvalueStudent(module.membership.measure.4, nSamples)


module.membership.measure.pvals.4[1:10,1:10]
#6
module.membership.measure.6 <- cor(module_eigengenes.6, norm.counts, use = 'p')
module.membership.measure.pvals.6 <- corPvalueStudent(module.membership.measure.6, nSamples)


module.membership.measure.pvals.6[1:10,1:10]
#8
module.membership.measure.8 <- cor(module_eigengenes.8, norm.counts, use = 'p')
module.membership.measure.pvals.8 <- corPvalueStudent(module.membership.measure.8, nSamples)


module.membership.measure.pvals.8[1:10,1:10]

###########################################
# Saving for the future and for TopGO analyses
###########################################

module.gene.mapping.4$geneID <- rownames(module.gene.mapping.4)
module.gene.mapping.6$geneID <- rownames(module.gene.mapping.6)
module.gene.mapping.8$geneID <- rownames(module.gene.mapping.8)

parts.4 <- strsplit(module.gene.mapping.4$geneID, "::")
parts.6 <- strsplit(module.gene.mapping.6$geneID, "::")
parts.8 <- strsplit(module.gene.mapping.8$geneID, "::")

module.gene.mapping.4$shortID <- sapply(
  parts.4,
  function(x) paste(x[1:2], collapse = "::")
)
module.gene.mapping.6$shortID <- sapply(
  parts.6,
  function(x) paste(x[1:2], collapse = "::")
)
module.gene.mapping.8$shortID <- sapply(
  parts.8,
  function(x) paste(x[1:2], collapse = "::")
)
write.csv(
  module.gene.mapping.4,
  "mydata/myresults/WGCNA_module_gene_mapping.4_w_shortIDs.csv",
  row.names = FALSE
)
write.csv(
  module.gene.mapping.6,
  "mydata/myresults/WGCNA_module_gene_mapping.6_w_shortIDs.csv",
  row.names = FALSE
)
write.csv(
  module.gene.mapping.8,
  "mydata/myresults/WGCNA_module_gene_mapping.8_w_shortIDs.csv",
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
  data.frame(shortID = yellow_shortIDs.4),
  "mydata/myresults/yellow_module_genes.4.csv",
  row.names = FALSE
)
#6
yellow_shortIDs.6 <- unique(
  module.gene.mapping.6$shortID[
    module.gene.mapping.6$`bwnet.6$colors` == "yellow"
  ]
)

length(yellow_shortIDs.6)

head(yellow_shortIDs.6)


write.csv(
  data.frame(shortID = yellow_shortIDs.6),
  "mydata/myresults/yellow_module_genes.6.csv",
  row.names = FALSE
)
#8
yellow_shortIDs.8 <- unique(
  module.gene.mapping.8$shortID[
    module.gene.mapping.8$`bwnet.8$colors` == "yellow"
  ]
)

length(yellow_shortIDs.8)

head(yellow_shortIDs.8)


write.csv(
  data.frame(shortID = yellow_shortIDs.8),
  "mydata/myresults/yellow_module_genes.8.csv",
  row.names = FALSE
)

turquoise_shortIDs.4 <- unique(
  module.gene.mapping.4$shortID[
    module.gene.mapping.4$`bwnet.4$colors` == "turquoise"
  ]
)

length(yellow_shortIDs.4)

head(yellow_shortIDs.4)


write.csv(
  data.frame(shortID = turquoise_shortIDs.4),
  "mydata/myresults/turquoise_module_genes.4.csv",
  row.names = FALSE
)
#6
turquoise_shortIDs.6 <- unique(
  module.gene.mapping.6$shortID[
    module.gene.mapping.6$`bwnet.6$colors` == "turquoise"
  ]
)

length(yellow_shortIDs.6)

head(yellow_shortIDs.6)


write.csv(
  data.frame(shortID = turquoise_shortIDs.6),
  "mydata/myresults/turquoise_module_genes.6.csv",
  row.names = FALSE
)
#8
turquoise_shortIDs.8 <- unique(
  module.gene.mapping.8$shortID[
    module.gene.mapping.8$`bwnet.8$colors` == "turquoise"
  ]
)

length(yellow_shortIDs.8)

head(yellow_shortIDs.8)


write.csv(
  data.frame(shortID = turquoise_shortIDs.8),
  "mydata/myresults/turquoise_module_genes.8.csv",
  row.names = FALSE
)





















turq4 <- rownames(module.gene.mapping.4)[module.gene.mapping.4$`bwnet.4$colors` == "turquoise"]
turq6 <- rownames(module.gene.mapping.6)[module.gene.mapping.6$`bwnet.6$colors` == "turquoise"]
turq8 <- rownames(module.gene.mapping.8)[module.gene.mapping.8$`bwnet.8$colors` == "turquoise"]


library(tidyverse)

sizes <- bind_rows(
  data.frame(power = 4, module = names(table(bwnet.4$colors)), n = as.integer(table(bwnet.4$colors))),
  data.frame(power = 6, module = names(table(bwnet.6$colors)), n = as.integer(table(bwnet.6$colors))),
  data.frame(power = 8, module = names(table(bwnet.8$colors)), n = as.integer(table(bwnet.8$colors)))
)

ggplot(sizes, aes(x = reorder(module, -n), y = n, fill = factor(power))) +
  geom_col(position = "dodge") +
  labs(x = "Module", y = "Number of genes", fill = "Soft power") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

table(bwnet.4$colors)
table(bwnet.6$colors)
table(bwnet.8$colors)

library(knitr)

module_counts <- data.frame(
  Module = c("turquoise", "blue", "brown", "grey", "yellow", "green", "red", "black", "pink", "magenta", "greenyellow", "purple", "tan", "salmon", "cyan"),
  Power4 = c(9221, 8226, 4676, 1878, 275, 234, 211, 200, 196, 143, 0, 0, 0, 0, 0),
  Power6 = c(9806, 7919, 3614, 1588, 1003, 367, 331, 233, 164, 144, 0, 91, 0, 0, 0),
  Power8 = c(7557, 6124, 4184, 3190, 1910, 676, 539, 325, 167, 162, 104, 114, 101, 61, 46)
)


##### GO Anaylsis Turquoise
module <- "turquoise"

# WGCNA module assignments
moduleColors <- bwnet.6$colors

# names(moduleColors) should be gene IDs
table(moduleColors)

geneList <- factor(
  as.integer(moduleColors == module)
)

names(geneList) <- sub(
  "::g\\.[0-9]+::m\\.[0-9]+$",
  "",
  names(bwnet.6$colors)
)

table(geneList)

##GO Analysis


module <- "turquoise"

# WGCNA module assignments
moduleColors.6 <- bwnet.6$colors

# names(moduleColors) should be gene IDs
table(moduleColors.6)

geneList <- factor(
  as.integer(moduleColors.6 == module)
)

names(geneList) <- sub(
  "::g\\.[0-9]+::m\\.[0-9]+$",
  "",
  names(bwnet.6$colors)
)

table(geneList)

##GO Analysis
library(topGO)

mappingFile <- "mydata/trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")

#Filter GO terms by size:(we dont want very small categories with small counts that 
#could be error, noise, or too specific, as well as super large ones (e.g. "biology" category))
allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]

GOdata <- new(
  "topGOdata",
  ontology = "BP",
  allGenes = geneList,
  geneSelectionFun = function(x) x == 1,
  annot = annFUN.gene2GO,
  gene2GO = geneID2GO.filtered
)

resultWeight <- runTest(
  GOdata,
  algorithm = "weight01",
  statistic = "fisher"
)

GOresults <- GenTable(
  GOdata,
  weightFisher = resultWeight,
  orderBy = "weightFisher",
  topNodes = 100
)

head(GOresults)
library(WGCNA)
library(GO.db)
vsd <- vst(dds)
datExpr <- t(assay(vsd))
MM <- cor(datExpr, bwnet.8$MEs, use = "p")

kME <- signedKME(
  datExpr,
  bwnet.8$MEs
)

hub.df <- data.frame(
  Gene = colnames(datExpr),
  Module = bwnet.8$colors,
  kME = kME$kMEpink
)

hub.df <- subset(
  hub.df,
  Module == "pink"
)

hub.df <- hub.df[
  order(-abs(hub.df$kME)),
]

head(hub.df, 5)

#Module       kME
#TRINITY_DN5073_c0_g1::TRINITY_DN5073_c0_g1_i1::g.20808::m.20808    pink 0.9317749
#TRINITY_DN3066_c0_g1::TRINITY_DN3066_c0_g1_i3::g.13912::m.13912    pink 0.9260290
#TRINITY_DN4365_c0_g1::TRINITY_DN4365_c0_g1_i3::g.18514::m.18514    pink 0.9200610
#TRINITY_DN4502_c0_g1::TRINITY_DN4502_c0_g1_i11::g.18894::m.18894   pink 0.9199457
#TRINITY_DN397_c36_g1::TRINITY_DN397_c36_g1_i1::g.2220::m.2220      pink 0.9185703
#grep("DN26848", names(geneID2GO), value = TRUE)

get_go <- function(id) {
  go <- geneID2GO[[id]]
  if (is.null(go)) return(data.frame(ID = id, GO.ID = NA, Term = NA, Ontology = NA))
  data.frame(ID = id, GO.ID = go, Term = Term(go), Ontology = Ontology(go))
}

top5 <- head(hub.df, 5)
top5$ID <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", rownames(top5))

top5$ID %in% names(geneID2GO)   # should all be TRUE now

top5_go <- do.call(rbind, lapply(top5$ID, get_go))
top5_go$kME <- top5$kME[match(top5_go$ID, top5$ID)]
top5_go revise this
# ============================================================
# 1. Rebuild datExpr / kME / hub.df (adjust bwnet.8 as needed)
# ============================================================
vsd     <- vst(dds)
datExpr <- t(assay(vsd))

kME <- WGCNA::signedKME(datExpr, bwnet.8$MEs)

hub.df <- data.frame(
  Gene   = colnames(datExpr),
  Module = bwnet.8$colors,
  kME    = kME$kMEpink,
  stringsAsFactors = FALSE
)

hub.df <- subset(hub.df, Module == "pink")
hub.df <- hub.df[order(-abs(hub.df$kME)), ]

top5 <- head(hub.df, 5)

# Extract the short ID (drop ::g.####::m.####)
top5$ID <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", top5$Gene)

cat("Top 5 hub genes (pink module):\n")
print(top5[, c("Gene", "ID", "kME")], row.names = FALSE)

cat("\nIDs found in geneID2GO:\n")
print(top5$ID %in% names(geneID2GO))

# ============================================================
# 2. Safe GO lookups (won't error on missing GO IDs)
# ============================================================
safe_term <- function(g) {
  t <- tryCatch(GO.db::Term(g),     error = function(e) NA_character_)
  if (length(t) == 0) NA_character_ else as.character(t[1])
}
safe_ont <- function(g) {
  o <- tryCatch(GO.db::Ontology(g), error = function(e) NA_character_)
  if (length(o) == 0) NA_character_ else as.character(o[1])
}

# ============================================================
# 3. get_go() — always returns the same columns, no rownames
# ============================================================
get_go <- function(id) {
  go <- geneID2GO[[id]]
  
  # Gene missing from the mapping → single NA row
  if (is.null(go) || length(go) == 0) {
    return(data.frame(
      ID       = id,
      GO.ID    = NA_character_,
      Term     = NA_character_,
      Ontology = NA_character_,
      row.names = NULL,
      stringsAsFactors = FALSE
    ))
  }
  
  # Gene present → one row per GO term
  data.frame(
    ID       = id,
    GO.ID    = as.character(go),
    Term     = vapply(go, safe_term, character(1)),
    Ontology = vapply(go, safe_ont,  character(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 4. Build top5_go — use bind_rows if available, else rbind
# ============================================================
go_rows <- lapply(top5$ID, get_go)

if (requireNamespace("dplyr", quietly = TRUE)) {
  top5_go <- dplyr::bind_rows(go_rows)
} else {
  # Ensure column names are identical across all pieces
  cols <- c("ID", "GO.ID", "Term", "Ontology")
  go_rows <- lapply(go_rows, function(x) x[, cols, drop = FALSE])
  rownames_all <- NULL
  top5_go <- do.call(rbind, c(go_rows, list(make.row.names = FALSE)))
  rownames(top5_go) <- NULL
}

# Attach kME per gene
top5_go$kME <- top5$kME[match(top5_go$ID, top5$ID)]

# ============================================================
# 5. Report
# ============================================================
cat("\nTop 5 hub gene GO annotations (pink module, P8):\n")
print(top5_go, row.names = FALSE)

# Optional: write to file
write.csv(top5_go,
          "mydata/myresults/top5_hub_GO_pink_P8.csv",
          row.names = FALSE)
cat("\nSaved: mydata/myresults/top5_hub_GO_pink_P8.csv\n")
#GO.ID                                        Term Annotated Significant Expected weightFisher
#1 GO:0006936                          muscle contraction       236         162    90.97      1.5e-26
#2 GO:0007527            adult somatic muscle development        85          74    32.77      1.2e-21
#3 GO:0071689               muscle thin filament assembly        81          70    31.22      3.9e-20
#4 GO:0043520 regulation of myosin II filament assembl...        74          65    28.52      1.3e-19
#5 GO:0030241 skeletal muscle myosin thick filament as...        86          72    33.15      5.1e-19
#6 GO:0060361                                      flight        83          68    31.99      4.9e-17

##Top Hub Genes:
# Calculate module membership
vsd <- vst(dds)
datExpr <- t(assay(vsd))
MM <- cor(datExpr, bwnet.8$MEs, use = "p")

kME <- signedKME(
  datExpr,
  bwnet.8$MEs
)

hub.df <- data.frame(
  Gene = colnames(datExpr),
  Module = bwnet.8$colors,
  kME = kME$kMEpink
)

hub.df <- subset(
  hub.df,
  Module == "pink"
)

hub.df <- hub.df[
  order(-abs(hub.df$kME)),
]

head(hub.df, 5)

#Module       kME
#TRINITY_DN5073_c0_g1::TRINITY_DN5073_c0_g1_i1::g.20808::m.20808    pink 0.9317749
#TRINITY_DN3066_c0_g1::TRINITY_DN3066_c0_g1_i3::g.13912::m.13912    pink 0.9260290
#TRINITY_DN4365_c0_g1::TRINITY_DN4365_c0_g1_i3::g.18514::m.18514    pink 0.9200610
#TRINITY_DN4502_c0_g1::TRINITY_DN4502_c0_g1_i11::g.18894::m.18894   pink 0.9199457
#TRINITY_DN397_c36_g1::TRINITY_DN397_c36_g1_i1::g.2220::m.2220      pink 0.9185703
#grep("DN26848", names(geneID2GO), value = TRUE)

get_go <- function(id) {
  go <- geneID2GO[[id]]
  if (is.null(go)) return(data.frame(ID = id, GO.ID = NA, Term = NA, Ontology = NA))
  data.frame(ID = id, GO.ID = go, Term = Term(go), Ontology = Ontology(go))
}

top5 <- head(hub.df, 5)
top5$ID <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", rownames(top5))

top5$ID %in% names(geneID2GO)   # should all be TRUE now

top5_go <- do.call(rbind, lapply(top5$ID, get_go))
top5_go$kME <- top5$kME[match(top5_go$ID, top5$ID)]
top5_go


#                                                           Term Ontology       kME
#1                                                               axon       CC 0.9668181
#2                                                 dense core granule       CC 0.9668181
#3                                                   exocytic vesicle       CC 0.9668181
#4                                                               <NA>     <NA> 0.9668181
#5                                                               <NA>     <NA> 0.9668181
#6                                             neuromuscular junction       CC 0.9668181
#7                                                    plasma membrane       CC 0.9668181
#8                                                   synaptic vesicle       CC 0.9668181
#9                                          synaptic vesicle membrane       CC 0.9668181
#10                                                   terminal bouton       CC 0.9668181

library(topGO)
library(dplyr)
library(tidyr)
library(ggplot2)

run_go <- function(net, power, module = "turquoise") {
  geneList <- factor(as.integer(net$colors == module))
  names(geneList) <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", names(net$colors))
  
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  res <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  # topNodes = all tested terms, so every term has a p-value at every power
  tab <- GenTable(GOdata, weightFisher = res, orderBy = "weightFisher",
                  topNodes = length(usedGO(GOdata)))
  
  tab$p     <- as.numeric(gsub("<", "", tab$weightFisher))
  tab$Power <- power
  tab
}

go4 <- run_go(bwnet.4, "4")
go6 <- run_go(bwnet.6, "6")
go8 <- run_go(bwnet.8, "8")

# Top 10 terms at each power
head(go4[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)
head(go6[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)
head(go8[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)










# 1. How many turquoise genes at each power, and how much do they overlap? ----
tq <- list(
  P4 = names(bwnet.4$colors)[bwnet.4$colors == "turquoise"],
  P6 = names(bwnet.6$colors)[bwnet.6$colors == "turquoise"],
  P8 = names(bwnet.8$colors)[bwnet.8$colors == "turquoise"]
)
sapply(tq, length)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
jaccard(tq$P4, tq$P6)
jaccard(tq$P4, tq$P8)
jaccard(tq$P6, tq$P8)

# 2. Overlap of the top-N enriched terms ---------------------------------------
N <- 50
top_terms <- list(
  P4 = head(go4$GO.ID, N),
  P6 = head(go6$GO.ID, N),
  P8 = head(go8$GO.ID, N)
)

length(intersect(top_terms$P4, top_terms$P6))
length(intersect(top_terms$P4, top_terms$P8))
length(intersect(top_terms$P6, top_terms$P8))
length(Reduce(intersect, top_terms))   # in the top N at all three powers

# 3. Side-by-side table of p-values --------------------------------------------
go_all <- bind_rows(go4, go6, go8)

go_wide <- go_all %>%
  dplyr::select(GO.ID, Term, Power, p) %>%
  pivot_wider(names_from = Power, values_from = p, names_prefix = "p_power") %>%
  mutate(best_p = pmin(p_power4, p_power6, p_power8, na.rm = TRUE)) %>%
  arrange(best_p)

head(go_wide, 20)

# Terms enriched at one power but not the others (p < 0.05 at only one power)
sig <- go_wide %>%
  dplyr::mutate(across(starts_with("p_power"), ~ .x < 0.05, .names = "sig_{.col}")) %>%
  dplyr::mutate(n_sig = sig_p_power4 + sig_p_power6 + sig_p_power8)

sig %>%
  dplyr::filter(n_sig == 1) %>%
  dplyr::select(GO.ID, Term, starts_with("p_power"))

# 4. Plot: top terms across powers ---------------------------------------------
top_union <- go_wide %>% slice_head(n = 15) %>% pull(GO.ID)

plot_df <- go_all %>%
  dplyr::filter(GO.ID %in% top_union) %>%
  dplyr::mutate(neglogp = -log10(pmax(p, 1e-30)),     # cap so tiny p's don't dominate
                Term = factor(Term, levels = rev(unique(go_wide$Term[go_wide$GO.ID %in% top_union]))))

p_go <- ggplot(plot_df, aes(x = neglogp, y = Term, fill = Power)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  theme_bw() +
  labs(title = "Top GO terms in the turquoise module across soft powers",
       x = "-log10(weightFisher p-value)", y = NULL, fill = "Soft power")







## Go Yellow



##GO Analysis


module <- "pink"

# WGCNA module assignments
moduleColors.8 <- bwnet.8$colors

# names(moduleColors) should be gene IDs
table(moduleColors.8)

geneListpink <- factor(
  as.integer(moduleColors.8 == module)
)

names(geneListpink) <- sub(
  "::g\\.[0-9]+::m\\.[0-9]+$",
  "",
  names(bwnet.8$colors)
)

table(geneListpink)

##GO Analysis
library(topGO)

mappingFile <- "mydata/trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")

#Filter GO terms by size:(we dont want very small categories with small counts that 
#could be error, noise, or too specific, as well as super large ones (e.g. "biology" category))
allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]
# Collapse duplicate gene entries: one list element per gene, all GO terms unioned
geneID2GO.filtered <- split(
  unlist(geneID2GO.filtered, use.names = FALSE),
  rep(names(geneID2GO.filtered), lengths(geneID2GO.filtered))
)

# Deduplicate GO terms within each gene
geneID2GO.filtered <- lapply(geneID2GO.filtered, unique)

# Drop empties
geneID2GO.filtered <- geneID2GO.filtered[lengths(geneID2GO.filtered) > 0]
GOdatapink <- new(
  "topGOdata",
  ontology = "BP",
  allGenes = geneListpink,
  geneSelectionFun = function(x) x == 1,
  annot = annFUN.gene2GO,
  gene2GO = geneID2GO.filtered
)

resultWeightpink <- runTest(
  GOdatapink,
  algorithm = "weight01",
  statistic = "fisher"
)

GOresultspink<- GenTable(
  GOdatapink,
  weightFisher = resultWeightpink,
  orderBy = "weightFisher",
  topNodes = 100
)

head(GOresultspink)

#GO.ID                                        Term Annotated Significant Expected
#1 GO:0010882 regulation of cardiac muscle contraction...        18           7     0.68
#2 GO:1901673      regulation of mitotic spindle assembly        10           7     0.38
#3 GO:0031291             Ran protein signal transduction         8           6     0.30
#4 GO:0042538              hyperosmotic salinity response         7           5     0.27
#5 GO:0030574                  collagen catabolic process        33           9     1.25
#6 GO:0045930 negative regulation of mitotic cell cycl...       152          10     5.78
#weightFisher
#1      2.6e-09
#2      1.0e-08
#3      6.9e-08
#4      1.4e-06
#5      2.3e-06
#6      7.3e-06

##Top Hub Genes:
# Calculate module membership
vsd <- vst(dds)
datExpr <- t(assay(vsd))
MM <- cor(datExpr, bwnet.8$MEs, use = "p")

kME <- signedKME(
  datExpr,
  bwnet.8$MEs
)

hub.df <- data.frame(
  Gene = colnames(datExpr),
  Module = bwnet.8$colors,
  kME = kME$kMEpink
)

hub.df.y <- subset(
  hub.df,
  Module == "pink"
)

hub.df.y <- hub.df[
  order(-abs(hub.df$kME)),
]

head(hub.df.y, 5)

#                                                                     Module       kME
#TRINITY_DN7934_c0_g2::TRINITY_DN7934_c0_g2_i2::g.28825::m.28825      yellow 0.9252903
#RINITY_DN194_c0_g1::TRINITY_DN194_c0_g1_i21::g.1129::m.1129      turquoise 0.9241266
#TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434    yellow 0.9140073#RINITY_DN5092_c0_g1::TRINITY_DN5092_c0_g1_i9::g.20841::m.20841      yellow 0.9067385
#TRINITY_DN5092_c0_g1::TRINITY_DN5092_c0_g1_i8::g.20839::m.20839      yellow 0.9065799
#> 
#grep("DN26848", names(geneID2GO), value = TRUE)

# ---- Safe lookups for GO Term / Ontology ----
safe_term <- function(g) {
  t <- tryCatch(GO.db::Term(g), error = function(e) NA_character_)
  if (length(t) == 0) NA_character_ else as.character(t[1])
}
safe_ont <- function(g) {
  o <- tryCatch(GO.db::Ontology(g), error = function(e) NA_character_)
  if (length(o) == 0) NA_character_ else as.character(o[1])
}

# ---- get_go() that always returns the same columns ----
get_go <- function(id) {
  go <- geneID2GO[[id]]
  if (is.null(go) || length(go) == 0) {
    return(data.frame(ID = id, GO.ID = NA_character_,
                      Term = NA_character_, Ontology = NA_character_,
                      stringsAsFactors = FALSE))
  }
  data.frame(
    ID       = id,
    GO.ID    = as.character(go),
    Term     = vapply(go, safe_term, character(1)),
    Ontology = vapply(go, safe_ont,  character(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ---- Build the top5 GO table ----
library(dplyr)
top5_go <- bind_rows(lapply(top5$ID, get_go))
top5_go$kME <- top5$kME[match(top5_go$ID, top5$ID)]
top5_go

#                                                           Term Ontology       kME
#1                                                               axon       CC 0.9668181
#2                                                 dense core granule       CC 0.9668181
#3                                                   exocytic vesicle       CC 0.9668181
#4                                                               <NA>     <NA> 0.9668181
#5                                                               <NA>     <NA> 0.9668181
#6                                             neuromuscular junction       CC 0.9668181
#7                                                    plasma membrane       CC 0.9668181
#8                                                   synaptic vesicle       CC 0.9668181
#9                                          synaptic vesicle membrane       CC 0.9668181
#10                                                   terminal bouton       CC 0.9668181

library(topGO)
library(dplyr)
library(tidyr)
library(ggplot2)

run_go_yellow <- function(net, power, module = "yellow") {
  geneList <- factor(as.integer(net$colors == module))
  names(geneList) <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", names(net$colors))
  
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  res <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  # topNodes = all tested terms, so every term has a p-value at every power
  tab <- GenTable(GOdata, weightFisher = res, orderBy = "weightFisher",
                  topNodes = length(usedGO(GOdata)))
  
  tab$p     <- as.numeric(gsub("<", "", tab$weightFisher))
  tab$Power <- power
  tab
}

go4y <- run_go_yellow(bwnet.4, "4")
go6y <- run_go_yellow(bwnet.6, "6")
go8y <- run_go_yellow(bwnet.8, "8")

# Top 10 terms at each power
head(go4y[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)
head(go6y[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)
head(go8y[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)



# 1. How many turquoise genes at each power, and how much do they overlap? ----
yl <- list(
  P4y = names(bwnet.4$colors)[bwnet.4$colors == "yellow"],
  P6y = names(bwnet.6$colors)[bwnet.6$colors == "yellow"],
  P8y = names(bwnet.8$colors)[bwnet.8$colors == "yellow"]
)
sapply(yl, length)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
jaccard(yl$P4y, yl$P6y)
jaccard(yl$P4y, yl$P8y)
jaccard(yl$P6y, yl$P8y)

# 2. Overlap of the top-N enriched terms ---------------------------------------
N <- 50
top_terms <- list(
  P4y = head(go4y$GO.ID, N),
  P6y = head(go6y$GO.ID, N),
  P8y = head(go8y$GO.ID, N)
)

length(intersect(top_terms$P4y, top_terms$P6y))
length(intersect(top_terms$P4y, top_terms$P8y))
length(intersect(top_terms$P6y, top_terms$P8y))
length(Reduce(intersect, top_terms))   # in the top N at all three powers
list(Reduce(intersect, top_terms))
# 3. Side-by-side table of p-values --------------------------------------------
go_all <- bind_rows(go4y, go6y, go8y)

go_wide <- go_all %>%
  dplyr::select(GO.ID, Term, Power, p) %>%
  pivot_wider(names_from = Power, values_from = p, names_prefix = "p_power") %>%
  mutate(best_p = pmin(p_power4, p_power6, p_power8, na.rm = TRUE)) %>%
  arrange(best_p)

head(go_wide, 20)

# Terms enriched at one power but not the others (p < 0.05 at only one power)
sig <- go_wide %>%
  dplyr::mutate(across(starts_with("p_power"), ~ .x < 0.05, .names = "sig_{.col}")) %>%
  dplyr::mutate(n_sig = sig_p_power4 + sig_p_power6 + sig_p_power8)

sig %>%
  dplyr::filter(n_sig == 1) %>%
  dplyr::select(GO.ID, Term, starts_with("p_power"))

# 4. Plot: top terms across powers ---------------------------------------------
top_union <- go_wide %>% slice_head(n = 15) %>% pull(GO.ID)

plot_df <- go_all %>%
  dplyr::filter(GO.ID %in% top_union) %>%
  dplyr::mutate(neglogp = -log10(pmax(p, 1e-30)),     # cap so tiny p's don't dominate
                Term = factor(Term, levels = rev(unique(go_wide$Term[go_wide$GO.ID %in% top_union]))))

p_go.y <- ggplot(plot_df, aes(x = neglogp, y = Term, fill = Power)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  theme_bw() +
  labs(title = "Top GO terms in the Yellow module",
       x = "-log10(weightFisher p-value)", y = NULL, fill = "Soft power")
print(p_go.y)
print(p_go)
## Fig 2

library(tidyverse)
library(patchwork)

# ---- A: all module sizes ----
grey_df <- data.frame(
  power = factor(c(4,6,8), levels = c(4,6,8)),
  grey  = c(sum(bwnet.4$colors == "grey"),
            sum(bwnet.6$colors == "grey"),
            sum(bwnet.8$colors == "grey")),
  total = c(length(bwnet.4$colors),
            length(bwnet.6$colors),
            length(bwnet.8$colors))
)
grey_df$pct <- 100 * grey_df$grey / grey_df$total

pA <- ggplot(grey_df, aes(power, pct, group = 1)) +
  geom_line(size = 1.2, color = "#555555") +
  geom_point(size = 5, color = "#555555") +
  geom_text(aes(label = paste0(round(pct,1), "%\n(", grey, ")")),
            vjust = -1.2, size = 4) +
  labs(x = "Soft power (β)",
       y = "Genes in Grey Module (%)",
       title = "A") +
  theme_bw(base_size = 12) +
  ylim(0, max(grey_df$pct) * 1.5)

pA

dom_df <- data.frame(
  power = factor(c(4,6,8), levels = c(4,6,8)),
  pct   = c(100 * max(table(bwnet.4$colors)) / length(bwnet.4$colors),
            100 * max(table(bwnet.6$colors)) / length(bwnet.6$colors),
            100 * max(table(bwnet.8$colors)) / length(bwnet.8$colors))
)

# ---- B: yellow + turquoise trajectories ----
focus <- data.frame(
  power  = factor(c(4,6,8, 4,6,8), levels = c(4,6,8)),
  module = rep(c("Turquoise","Yellow"), each = 3),
  n      = c(9221, 9806, 7557, 275, 1003, 1910)
)
pB <- ggplot(focus, aes(power, n, color = module, group = module)) +
  geom_line(size = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = c("Turquoise" = "#40E0D0",
                                "Yellow"    = "#E6B800")) +
  labs(x = "Soft power (β)", y = "Number of genes", title = "B") +
  theme_bw()

# ---- C: Jaccard heatmap ----
jac_df <- data.frame(
  module  = rep(c("Turquoise","Yellow"), each = 3),
  pair    = rep(c("4 vs 6","4 vs 8","6 vs 8"), times = 2),
  jaccard = c(
    jaccard(tq$P4, tq$P6), jaccard(tq$P4, tq$P8), jaccard(tq$P6, tq$P8),
    jaccard(yl$P4y, yl$P6y), jaccard(yl$P4y, yl$P8y), jaccard(yl$P6y, yl$P8y)
  )
)
pC <- ggplot(jac_df, aes(pair, module, fill = jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(jaccard, 2)), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "#CC3333", limits = c(0,1)) +
  labs(x = "Power comparison", y = NULL, title = "C", fill = "Jaccard") +
  theme_minimal()

# Combine
(pA / (pB | pC)) +
  plot_annotation(
    title = ""
  )

library(WGCNA)
library(CorLevelPlot)
library(dplyr)
library(tidyr)

# ---- Traits ----
traits <- data.frame(
  Warming       = as.numeric(sample_metadata$treatment %in% c("OW", "OWA")),
  Acidification = as.numeric(sample_metadata$treatment %in% c("OA", "OWA")),
  WarmingAcidification   = as.numeric(sample_metadata$treatment == "OWA"),
  row.names     = rownames(sample_metadata)
)

# ---- P4 ----
{
  MEs <- orderMEs(bwnet.4$MEs)
  heatmap.data <- merge(MEs, traits, by = "row.names") %>%
    column_to_rownames("Row.names")
  
  pdf("mydata/myresults/module_trait_P4.pdf", width = 8, height = 7)
  print(CorLevelPlot(heatmap.data,
                     x = colnames(traits), y = colnames(MEs),
                     col = c("blue1", "skyblue", "white", "pink", "red"),
                     main = "Module-Trait Relationships (power P4)"))
  dev.off()
  cat("Wrote P4\n")
}

# ---- P6 ----
{
  MEs <- orderMEs(bwnet.6$MEs)
  heatmap.data <- merge(MEs, traits, by = "row.names") %>%
    column_to_rownames("Row.names")
  
  pdf("mydata/myresults/module_trait_P6.pdf", width = 8, height = 7)
  print(CorLevelPlot(heatmap.data,
                     x = colnames(traits), y = colnames(MEs),
                     col = c("blue1", "skyblue", "white", "pink", "red"),
                     main = "Module-Trait Relationships (power P6)"))
  dev.off()
  cat("Wrote P6\n")
}

# ---- P8 ----
{
  MEs <- orderMEs(bwnet.8$MEs)
  heatmap.data <- merge(MEs, traits, by = "row.names") %>%
    column_to_rownames("Row.names")
  
  pdf("mydata/myresults/module_trait_P8.pdf", width = 8, height = 7)
  print(CorLevelPlot(heatmap.data,
                     x = colnames(traits), y = colnames(MEs),
                     col = c("blue1", "skyblue", "white", "pink", "red"),
                     main = "Module-Trait Relationships (power P8)"))
  dev.off()
  cat("Wrote P8\n")
}

# =============================================================
# kME summary per module, at each soft power
# Base R only — no dplyr, no n() conflicts
# =============================================================

library(WGCNA)

# ---- 0. Make sure datExpr exists ----
if (!exists("datExpr")) {
  vsd <- vst(dds)
  datExpr <- t(assay(vsd))
}

# ---- 1. Define the networks ----
nets <- list(P4 = bwnet.4, P6 = bwnet.6, P8 = bwnet.8)

# ---- 2. Function: kME stats for one network ----
compute_kME_stats <- function(net, power_label, expr) {
  cat("\n>>> Computing kME for", power_label, "...\n")
  
  # signedKME returns a genes x modules matrix
  kme <- WGCNA::signedKME(expr, net$MEs)
  
  # Self-kME: for each gene, the kME of its OWN module
  module_vec  <- net$colors
  self_colname <- paste0("kME", module_vec)
  self_kME <- kme[cbind(seq_len(nrow(kme)),
                        match(self_colname, colnames(kme)))]
  
  # Build per-gene data frame
  per_gene <- data.frame(
    gene     = rownames(kme),
    module   = module_vec,
    kME_self = as.numeric(self_kME),
    kme,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  
  # ---- Per-module summary via base R ----
  modules <- sort(unique(module_vec))
  summary_rows <- lapply(modules, function(m) {
    x <- per_gene$kME_self[per_gene$module == m]
    data.frame(
      module      = m,
      n           = length(x),
      mean_kME    = round(mean(x),   3),
      median_kME  = round(median(x), 3),
      sd_kME      = round(sd(x),     3),
      q25_kME     = round(quantile(x, 0.25), 3),
      q75_kME     = round(quantile(x, 0.75), 3),
      frac_gt_0.7 = round(mean(x > 0.7), 2),
      power       = power_label,
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(-summary_df$median_kME), ]
  rownames(summary_df) <- NULL
  
  list(per_gene = per_gene, summary = summary_df)
}

# ---- 3. Run for all three powers ----
results <- list()
for (nm in names(nets)) {
  results[[nm]] <- compute_kME_stats(nets[[nm]], nm, datExpr)
}

# ---- 4. Combine summaries ----
kme_summary_all <- do.call(rbind, lapply(results, `[[`, "summary"))
rownames(kme_summary_all) <- NULL

# Save
write.csv(kme_summary_all,
          "mydata/myresults/kME_summary_all_powers.csv",
          row.names = FALSE)

# Save per-gene tables
for (nm in names(results)) {
  write.csv(results[[nm]]$per_gene,
            paste0("mydata/myresults/kME_pergene_", nm, ".csv"),
            row.names = FALSE)
}

# ---- 5. Print per-power summary ----
for (nm in names(results)) {
  cat("\n===========================================\n")
  cat("Median kME —", nm, "\n")
  cat("===========================================\n")
  print(results[[nm]]$summary, row.names = FALSE)
}

# ---- 6. Side-by-side across powers ----
modules_all <- sort(unique(kme_summary_all$module))
powers_all  <- sort(unique(kme_summary_all$power))

kme_wide <- data.frame(module = modules_all, stringsAsFactors = FALSE)
for (p in powers_all) {
  sub <- kme_summary_all[kme_summary_all$power == p, c("module", "median_kME")]
  kme_wide[[paste0("median_kME_", p)]] <-
    sub$median_kME[match(kme_wide$module, sub$module)]
}

# Sort by P8 (last power) descending, NAs last
kme_wide <- kme_wide[order(-kme_wide[[paste0("median_kME_", powers_all[length(powers_all)]])],
                           na.last = TRUE), ]
rownames(kme_wide) <- NULL

cat("\n===========================================\n")
cat("Median kME across powers\n")
cat("===========================================\n")
print(kme_wide, row.names = FALSE)

write.csv(kme_wide,
          "mydata/myresults/kME_median_wide.csv",
          row.names = FALSE)

cat("\nFiles written:\n")
cat("  mydata/myresults/kME_summary_all_powers.csv\n")
cat("  mydata/myresults/kME_median_wide.csv\n")
cat("  mydata/myresults/kME_pergene_P4.csv\n")
cat("  mydata/myresults/kME_pergene_P6.csv\n")
cat("  mydata/myresults/kME_pergene_P8.csv\n")

 install.packages("gt")   # if not installed
 library(gt)
 library(dplyr)
 library(scales)
 
 # ---- Your data ----
 kme_wide <- data.frame(
   module        = c("black","blue","brown","cyan","green","greenyellow",
                     "grey","magenta","pink","purple","red","salmon","tan",
                     "turquoise","yellow"),
   median_kME_P4 = c(0.596,0.568,0.519,NA,0.750,NA,0.038,0.747,0.773,NA,
                     0.657,NA,NA,0.576,0.514),
   median_kME_P6 = c(0.640,0.579,0.582,NA,-0.353,NA,0.055,0.599,0.776,0.774,
                     0.649,NA,NA,0.539,0.609),
   median_kME_P8 = c(0.755,0.666,0.571,0.817,0.725,0.804,0.023,0.498,0.768,
                     0.739,0.650,0.803,0.813,0.637,0.446),
   stringsAsFactors = FALSE
 )
 
 # ---- Build the gt table ----
 tbl <- kme_wide %>%
   gt() %>%
   tab_header(
     title    = "Median kME Across Soft Powers",
     subtitle = ""
   ) %>%
   cols_label(
     module        = "Module",
     median_kME_P4 = "Power 4",
     median_kME_P6 = "Power 6",
     median_kME_P8 = "Power 8"
   ) %>%
   fmt_number(
     columns  = starts_with("median_kME"),
     decimals = 2
   ) %>%
   sub_missing(
     columns      = starts_with("median_kME"),
     missing_text = "—"
   ) %>%
   data_color(
     columns = starts_with("median_kME"),
     fn = scales::col_numeric(
       palette  = c("steelblue", "white", "firebrick"),
       domain   = c(-0.5, 1),          # ← widened to include negative green
       na.color = "grey90"
     )
   ) %>%
   cols_align(align = "center", columns = everything()) %>%
   tab_options(
     table.font.size         = 12,
     heading.title.font.size = 14,
     heading.align           = "left"
   )
 
 print(tbl)
 
 gtsave(tbl, "mydata/myresults/kME_median_table.html")
 