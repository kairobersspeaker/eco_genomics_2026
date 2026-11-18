library(tidyverse)
OUTDIR="/gpfs1/home/k/i/kirobert/Projects/eco_genomics_2026/population_genomics/myresults/HW"
# Read in population data
pops <- readLines("baypass/Lottia.baypass.pops.txt")
# define number of poulations object
npop <- length(pops)
# Read in environmental data per population
env <- read_csv("pop_level_biooracle.csv")
# order populations by lat (south to north)
pop_order<-env$pop[order(env$lat)]
# Read in omega matrix
omega <- as.matrix(read.table("BayPassHW/Lottia_bp_2_mat_omega.out"))
# Define column/row values of omega matrix by appropriate population 
dimnames(omega) <- list(pops, pops)
#scales omega vales to be on comparable scale
omega_cor <- cov2cor(omega)
# prints human readable matrix for validation
round(omega_cor, 2)
#Build heirarchical clustering dendrogram
omega_tree <- hclust(as.dist(1 - omega_cor), method = "average")
## Generate and save omega tree
png(paste0(OUTDIR,"/baypass_omega_tree.png"), width = 1800, height = 1200, res = 300)
plot(omega_tree, main = "Populations clustered by BayPass Omega", xlab = "", sub = "")
dev.off()
#Reshape omega matrix for usage in heatmap
omega_df <- as.data.frame(as.table(omega_cor)) %>%
  rename(pop1 = Var1, pop2 = Var2, correlation = Freq) %>%
  mutate(pop1 = factor(pop1, levels = pop_order),
         pop2 = factor(pop2, levels = pop_order))
###plot and save heatmap
ggplot(omega_df, aes(x = pop1, y = pop2, fill = correlation)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "cadetblue", mid = "white", high = "red4",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL,
       title = "BayPass Omega as correlations (south to north)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(OUTDIR,"/baypass_omega.png"), width = 7, height = 6)
# read in BayPass snps.txt metadata file
snps <- read_tsv("baypass/Lottia.baypass.snps.txt")
# assigns numeric row marker column 
snps$MRK<-rownames(snps)
snps$MRK<-as.numeric(snps$MRK)
# Defines population number
npop <- 19
# Read in BayPass XtX output file
xtx <- read.table("BayPassHW/Lottia_bp_2_summary_pi_xtx.out", header = TRUE, check.names = FALSE) %>%
  left_join(snps, by = "MRK") %>%
  mutate(lp_high   = pchisq(XtXst, npop, lower.tail = FALSE, log.p = TRUE),   #
         lp_low    = pchisq(XtXst, npop, lower.tail = TRUE,  log.p = TRUE),   
         direction = if_else(lp_high < lp_low, "high", "low"),
         log10_p   = -(pmin(lp_high, lp_low) + log(2)) / log(10),             
         pval      = 10^(-log10_p),
         qval      = p.adjust(pval, method = "BH"))
#plot XtX p-val for data validation
ggplot(xtx, aes(x = pval)) +
  geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "grey60", colour = "white") +
  labs(x = "p-value", y = "no. of SNPs") +
  theme_classic()
ggsave(paste0(OUTDIR,"/baypass_xtx_pvalues.png"), width = 5, height = 4)
##Adding genome position to SNPs
add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}
#apply position to XtX object
xtx <- xtx %>%
  add_genome_pos()    
# Define threshold :top 0.1% (smallest q) 
threshold <- quantile(xtx$qval, 0.001)   
#create ouliers object as a function of threshold and XtX
xtx_outliers <- xtx %>%
  filter(qval < threshold) %>%
  arrange(desc(qval)) %>%
  select(MRK, CHROM, POS, M_XtX, XtXst, pval, qval)
#observe and save outliers
head(xtx_outliers, 10)
write_tsv(xtx_outliers, "BayPassHW/Lottia_bp_xtx_outliers.txt")
#Load and Manipulate covariate output file for anaylisis 
bf <- read.table("BayPassHW/Lottia_bp_2_summary_betai_reg.out", header = TRUE, check.names = FALSE) %>%
  filter(COVARIABLE %in% c(2, 4)) %>%
  rename(BF_dB = `BF(dB)`) %>%
  mutate(covariate = ifelse(COVARIABLE == 2, "O2mean", "Chl_mean")) %>%
  left_join(snps, by = "MRK") %>%
  add_genome_pos()
#Counts SNPs with decisive (x>20) bayes factor
bf %>%
  group_by(covariate) %>%
  summarise(n_decisive = sum(BF_dB > 20))
# A tibble: 2 × 2
#covariate n_decisive
#<chr>          <int>
#  1 Chl_mean          76
#2 O2mean           333

# 3. Plot and save the Bayes Factor Scan
ggplot(bf, aes(x = genome_pos / 1e6, y = BF_dB)) +
  geom_point(aes(colour = shade), size = 0.6) +
  geom_point(data = filter(bf, BF_dB > 20), colour = "red", size = 1) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  facet_wrap(~ covariate, ncol = 1) +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Genome pos (Mb)",
       y = "Bayes factor (dB)") +
  theme_classic()
ggsave(paste0(OUTDIR,"/baypass_bf_scan_HW.png"), width = 10, height = 6)
# read in geno file
geno <- as.matrix(read.table("baypass/Lottia.baypass.geno"))
plot_top_snp <- function(target_row, top_5_snps,
                         tag,          # short label for title/files: "Chl" or "O2"
                         env_var,      # env column name: "chl_mean" or "o2_mean"
                         xlab) {       # full x-axis label incl. units
  
  # Pick the Nth-ranked SNP and pull its genotype row
  top_snp <- top_5_snps[target_row, ]
  counts  <- geno[top_snp$MRK, ]
  
  # BayPass .geno rows alternate ref, alt across populations 
  ref <- counts[seq(1, length(counts), by = 2)]
  alt <- counts[seq(2, length(counts), by = 2)]
  # Per population alt allele frequency joined with environmental values
  top_freq <- tibble(pop = pops, alt_freq = alt / (ref + alt)) %>%
    left_join(env, by = "pop")
   # Scatter of allele frequency vs. environment, points labelled by population
  p <- ggplot(top_freq, aes(x = .data[[env_var]], y = alt_freq, label = pop)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    ylim(0, 1) +
    # Title: short tag + rank + chromosome:position
    labs(x = xlab, y = "Allele Frequency",
         title = paste0(tag, " ", target_row, ": ",
                        top_snp$CHROM, ":", top_snp$POS)) +
    theme_classic()
  
  # Save per-rank PNG (filename includes tag + rank), then return the plot
  ggsave(paste0(OUTDIR, "/baypass_top_", tag, "_snp_rank", target_row, ".png"),
         p, width = 5, height = 4)
  p
}

# --- Chlorophyll: top 5 SNPs, then loop over ranks 1-5 ---
top_5_chl <- bf %>% filter(covariate == "Chl_mean") %>%
  slice_max(BF_dB, n = 5, with_ties = FALSE)   # top 5 by BF, no ties

chl_plots <- purrr::map(1:5, ~ plot_top_snp(   # run function for each rank
  .x, top_5_chl,
  tag     = "Chl",
  env_var = "chl_mean",
  xlab    = "Mean Chlorophyll Concentration (mg/m^3)"
))

# --- Oxygen: same steps for the O2 covariate ---
top_5_o2 <- bf %>% filter(covariate == "O2mean") %>%
  slice_max(BF_dB, n = 5, with_ties = FALSE)

o2_plots <- purrr::map(1:5, ~ plot_top_snp(
  .x, top_5_o2,
  tag     = "O2",
  env_var = "o2_mean",
  xlab    = "Mean O2 Concentration (mmol/m^3)"
))
####annotation####
# 1. Load the SnpEff gene effects file
geneanno <- read.table("outlier_gene_effects.tsv")
names(geneanno) <- c("CHROM", "POS", "GeneID", "Annotation")

# 2. Get your outliers (BF > 20) for O2 and Chl
my_outliers <- bf %>%
  filter(BF_dB > 20, covariate %in% c("O2mean", "Chl_mean")) %>%
  dplyr::select(CHROM, POS, covariate, BF_dB)
# 3. Match them to genes
my_outliers_annotated <- my_outliers %>%
  left_join(geneanno, by = c("CHROM", "POS"), relationship = "many-to-many") %>%
  filter(!is.na(Annotation)) # Keep only those inside/next to genes

# 4. See the top 5 most strongly associated genes for each variable
top_genes <- my_outliers_annotated %>%
  group_by(covariate) %>%
  arrange(desc(BF_dB)) %>%
  slice_head(n = 5)

print(top_genes)

# 5. Plot the types of mutations (matches Section 4.2 in your doc)
ggplot(my_outliers_annotated, aes(x = Annotation, fill = Annotation)) + 
  geom_bar(stat = "count") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~covariate) +
  labs(title = "Types of mutations in environmental outlier SNPs", x = "Effect", y = "Count")

ggsave(paste0(OUTDIR, "/baypass_snp_effects.png"), width = 8, height = 6)




####GO####
library(topGO)

# 1. Get the candidate genes for O2 ONLY
universe <- names(gene2GO)
candidates_o2 <- my_outliers %>%
  filter(covariate == "O2mean") %>%
  left_join(snp2gene, by = c("CHROM", "POS"), relationship = "many-to-many") %>%
  pull(GeneId) %>%
  unique()

cat("Found", length(candidates_o2), "candidate genes for O2.\n")

# 2. Set up the gene list (1 = candidate, 0 = background)
geneList <- factor(as.integer(universe %in% candidates_o2), levels = c(0, 1))
names(geneList) <- universe

# 3. Build the topGO object
GOdata_o2 <- new("topGOdata",
                 ontology = "BP", # Biological Process
                 allGenes = geneList,
                 annot = annFUN.gene2GO,
                 gene2GO = gene2GO,
                 nodeSize = 10)

# 4. Run the Fisher test
result_o2 <- runTest(GOdata_o2, algorithm = "weight01", statistic = "fisher")
res_o2 <- GenTable(GOdata_o2, weight01 = result_o2, orderBy = "weight01",
                   topNodes = 20, numChar = 1000)

# Clean up the p-value column
res_o2$p <- as.numeric(sub("< ", "", res_o2$weight01)) 

# Print the top 10 results to the console for your homework notes
print(head(res_o2[, c("Term", "Annotated", "Significant", "Expected", "p")], 10))
#Term Annotated Significant
#1  homophilic cell adhesion via plasma membrane adhesion molecules        99           5
#2                            post-Golgi vesicle-mediated transport        32           2
#3                                               protein processing        48           2
#4                                       response to light stimulus        10           1
#5                                                  DNA replication       250           4
#6                         cortical actin cytoskeleton organization        12           1
#7                                                 circadian rhythm        12           1
#8         double-strand break repair via break-induced replication        12           1
#9                         homologous chromosome pairing at meiosis        14           1
#10                                        organic cation transport        14           1
#Expected        p
#1      0.43 0.000067
#2      0.14 0.008500
#3      0.21 0.018600
#4      0.04 0.043000
#5      1.09 0.043700
#6      0.05 0.051300
#7      0.05 0.051300
#8      0.05 0.051300
#9      0.06 0.059600
#10     0.06 0.059600
#Now For Chl
candidates_chl <- my_outliers %>%
  filter(covariate == "Chl_mean") %>%
  left_join(snp2gene, by = c("CHROM", "POS"), relationship = "many-to-many") %>%
  pull(GeneId) %>%
  unique()

cat("Found", length(candidates_chl), "candidate genes for chl.\n")

# 2. Set up the gene list (1 = candidate, 0 = background)
geneList <- factor(as.integer(universe %in% candidates_chl), levels = c(0, 1))
names(geneList) <- universe

# 3. Build the topGO object
GOdata_chl <- new("topGOdata",
                  ontology = "BP", # Biological Process
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO,
                  nodeSize = 10)

# 4. Run the Fisher test
result_chl <- runTest(GOdata_chl, algorithm = "weight01", statistic = "fisher")
res_chl <- GenTable(GOdata_chl, weight01 = result_chl, orderBy = "weight01",
                    topNodes = 20, numChar = 1000)

# Clean up the p-value column
res_chl$p <- as.numeric(sub("< ", "", res_chl$weight01)) 

# Print the top 10 results to the console for your homework notes
print(head(res_chl[, c("Term", "Annotated", "Significant", "Expected", "p")], 10))
Term Annotated Significant
#1                                SMAD protein signal transduction        15           1
#2                                           BMP signaling pathway        34           1
#3                          negative regulation of gene expression       122           1
#4                              anatomical structure morphogenesis       149           1
#5                                            cell differentiation       317           1
#6                regulation of transcription by RNA polymerase II       756           1
#7                                                       signaling      1594           1
#8                                anatomical structure development       432           1
#9                          regulation of RNA biosynthetic process      1052           1
#10 regulation of nucleobase-containing compound metabolic process      1154           1
#Expected      p
#1      0.00 0.0031
#2      0.01 0.0071
#3      0.03 0.0253
#4      0.03 0.0308
#5      0.07 0.0650
#6      0.16 0.1514
#7      0.33 1.0000
#8      0.09 1.0000
#9      0.22 1.0000
#10     0.24 1.0000

