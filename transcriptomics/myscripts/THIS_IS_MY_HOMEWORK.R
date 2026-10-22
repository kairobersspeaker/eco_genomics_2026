## Set your working directory

# Load the packages
library(WGCNA)
options(stringsAsFactors = FALSE)

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)
library(pheatmap)
library(topGO)
library(dplyr)
library(tidyr)
library(patchwork)
library(GO.db)

# ============================================================
# 1. Import counts matrix and metadata, filter using DESeq2
# ============================================================

countsTable <- read.table("mydata/salmon.isoform.counts.matrix.filteredAssembly",
                          header = TRUE, row.names = 1)
countsTableRound <- round(countsTable)

sample_metadata <- read.table(file = "mydata/Ahud_trait_data.txt",
                              header = TRUE, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound,
                              colData   = sample_metadata,
                              design    = ~ 1)

# Filter: keep genes with >= 15 counts in >= 28 samples
dds <- dds[rowSums(counts(dds) >= 15) >= 28, ]
nrow(dds)

# Run DESeq model
dds <- DESeq(dds)

# ============================================================
# 2. QC - outlier detection
# ============================================================

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# Hierarchical clustering of samples
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree)

# PCA
pca <- prcomp(t(countsTable))
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, digits = 2)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# ============================================================
# 3. Normalization
# ============================================================

colData <- row.names(sample_metadata)

all(rownames(colData) %in% colnames(countsTableRound))
all(rownames(colData) == colnames(countsTableRound))

# Variance stabilization
dds_norm <- vst(dds)
norm.counts <- assay(dds_norm) %>% t()

# ============================================================
# 4. Network Construction
# ============================================================

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

# Convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

# Load pre-computed networks (blockwiseModules is memory-intensive)
bwnet.4 <- readRDS("mydata/bwnet.4.rds")
bwnet.6 <- readRDS("mydata/bwnet.6.rds")
bwnet.8 <- readRDS("mydata/bwnet.8.rds")

# ============================================================
# 5. Module Eigengenes
# ============================================================

module_eigengenes.4 <- bwnet.4$MEs
module_eigengenes.6 <- bwnet.6$MEs
module_eigengenes.8 <- bwnet.8$MEs

table(bwnet.4$colors)
table(bwnet.6$colors)
table(bwnet.8$colors)

# ============================================================
# 6. Dendrograms
# ============================================================

# P4
plotDendroAndColors(bwnet.4$dendrograms[[1]],
                    cbind(bwnet.4$unmergedColors, bwnet.4$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# P6
plotDendroAndColors(bwnet.6$dendrograms[[1]],
                    cbind(bwnet.6$unmergedColors, bwnet.6$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# P8
plotDendroAndColors(bwnet.8$dendrograms[[1]],
                    cbind(bwnet.8$unmergedColors, bwnet.8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# ============================================================
# 7A. Module–Trait Relationships (raw metadata columns)
# ============================================================

library(CorLevelPlot)
library(dplyr)
library(tidyr)

# Raw trait columns from sample_metadata (same as original script)
traits <- sample_metadata[, c(5, 8, 11, 14, 17)]

nSamples <- nrow(norm.counts)
nGenes   <- ncol(norm.counts)

trait_names <- colnames(traits)

# --- P4 ------------------------------------------------------------
MEs.4          <- orderMEs(bwnet.4$MEs)
heatmap.data.4 <- merge(MEs.4, traits, by = "row.names") %>%
  column_to_rownames("Row.names")

CorLevelPlot(heatmap.data.4,
             x = trait_names,
             y = colnames(MEs.4),
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = "Module–Trait Relationships (power 4)")

# --- P6 ------------------------------------------------------------
MEs.6          <- orderMEs(bwnet.6$MEs)
heatmap.data.6 <- merge(MEs.6, traits, by = "row.names") %>%
  column_to_rownames("Row.names")

CorLevelPlot(heatmap.data.6,
             x = trait_names,
             y = colnames(MEs.6),
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = "Module–Trait Relationships (power 6)")

# --- P8 ------------------------------------------------------------
MEs.8          <- orderMEs(bwnet.8$MEs)
heatmap.data.8 <- merge(MEs.8, traits, by = "row.names") %>%
  column_to_rownames("Row.names")

CorLevelPlot(heatmap.data.8,
             x = trait_names,
             y = colnames(MEs.8),
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = "Module–Trait Relationships (power 8)")

# ============================================================
# 7. Module–Treatment Correlations (heatmaps at each power)
# ============================================================

library(CorLevelPlot)
library(dplyr)
library(tidyr)

# --- Build binary treatment traits from sample_metadata$treatment ---
# Adjust the string values to match your metadata
# (e.g. "Ambient", "OA", "OW", "OWA" — check with table(sample_metadata$treatment))
traits <- data.frame(
  Warming             = as.numeric(sample_metadata$treatment %in% c("OW", "OWA")),
  Acidification       = as.numeric(sample_metadata$treatment %in% c("OA", "OWA")),
  `Warming/Acidification`    = as.numeric(sample_metadata$treatment == "OWA"),
  row.names           = rownames(sample_metadata)
)

# --- P4 ------------------------------------------------------------
MEs.4        <- orderMEs(bwnet.4$MEs)
heatmap.data.4 <- merge(MEs.4, traits, by = "row.names") %>%
  column_to_rownames(var = "Row.names")

CorLevelPlot(heatmap.data.4,
             x = colnames(traits),
             y = colnames(MEs.4),
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = "Module–Treatment Relationships (power 4)")

# --- P6 ------------------------------------------------------------
MEs.6        <- orderMEs(bwnet.6$MEs)
heatmap.data.6 <- merge(MEs.6, traits, by = "row.names") %>%
  column_to_rownames(var = "Row.names")

CorLevelPlot(heatmap.data.6,
             x = colnames(traits),
             y = colnames(MEs.6),
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = "Module–Treatment Relationships (power 6)")

# --- P8 ------------------------------------------------------------
MEs.8        <- orderMEs(bwnet.8$MEs)
heatmap.data.8 <- merge(MEs.8, traits, by = "row.names") %>%
  column_to_rownames(var = "Row.names")

CorLevelPlot(heatmap.data.8,
             x = colnames(traits),
             y = colnames(MEs.8),
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = "Module–Treatment Relationships (power 8)")
# ============================================================
# 8. Module Gene Mappings
# ============================================================

module.gene.mapping.4 <- as.data.frame(bwnet.4$colors)
module.gene.mapping.6 <- as.data.frame(bwnet.6$colors)
module.gene.mapping.8 <- as.data.frame(bwnet.8$colors)

# ============================================================
# 9. Figure 2: Module sizes / power comparison / Jaccard
# ============================================================

# ---- A: Grey module percentage ----
grey_df <- data.frame(
  power = factor(c(4, 6, 8), levels = c(4, 6, 8)),
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
  geom_text(aes(label = paste0(round(pct, 1), "%\n(", grey, ")")),
            vjust = -1.2, size = 4) +
  labs(x = "Soft power (β)",
       y = "Genes in Grey Module (%)",
       title = "A") +
  theme_bw(base_size = 12) +
  ylim(0, max(grey_df$pct) * 1.5)

# ---- B: Yellow + turquoise trajectories ----
focus <- data.frame(
  power  = factor(c(4, 6, 8, 4, 6, 8), levels = c(4, 6, 8)),
  module = rep(c("Turquoise", "Yellow"), each = 3),
  n      = c(9221, 9806, 7557, 275, 1003, 1910)
)

pB <- ggplot(focus, aes(power, n, color = module, group = module)) +
  geom_line(size = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Turquoise" = "#40E0D0",
                                "Yellow"    = "#E6B800")) +
  labs(x = "Soft power (β)", y = "Number of genes", title = "B") +
  theme_bw()

# ---- C: Jaccard heatmap ----
tq <- list(
  P4 = names(bwnet.4$colors)[bwnet.4$colors == "turquoise"],
  P6 = names(bwnet.6$colors)[bwnet.6$colors == "turquoise"],
  P8 = names(bwnet.8$colors)[bwnet.8$colors == "turquoise"]
)
yl <- list(
  P4y = names(bwnet.4$colors)[bwnet.4$colors == "yellow"],
  P6y = names(bwnet.6$colors)[bwnet.6$colors == "yellow"],
  P8y = names(bwnet.8$colors)[bwnet.8$colors == "yellow"]
)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

jac_df <- data.frame(
  module  = rep(c("Turquoise", "Yellow"), each = 3),
  pair    = rep(c("4 vs 6", "4 vs 8", "6 vs 8"), times = 2),
  jaccard = c(
    jaccard(tq$P4, tq$P6), jaccard(tq$P4, tq$P8), jaccard(tq$P6, tq$P8),
    jaccard(yl$P4y, yl$P6y), jaccard(yl$P4y, yl$P8y), jaccard(yl$P6y, yl$P8y)
  )
)

pC <- ggplot(jac_df, aes(pair, module, fill = jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(jaccard, 2)), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "#CC3333", limits = c(0, 1)) +
  labs(x = "Power comparison", y = NULL, title = "C", fill = "Jaccard") +
  theme_minimal()

# Combine Figure 2
(pA / (pB | pC)) +
  plot_annotation(title = "")

# ============================================================
# 10. Barplot: Module sizes across powers
# ============================================================

sizes <- bind_rows(
  data.frame(power = 4, module = names(table(bwnet.4$colors)),
             n = as.integer(table(bwnet.4$colors))),
  data.frame(power = 6, module = names(table(bwnet.6$colors)),
             n = as.integer(table(bwnet.6$colors))),
  data.frame(power = 8, module = names(table(bwnet.8$colors)),
             n = as.integer(table(bwnet.8$colors)))
)

ggplot(sizes, aes(x = reorder(module, -n), y = n, fill = factor(power))) +
  geom_col(position = "dodge") +
  labs(x = "Module", y = "Number of genes", fill = "Soft power") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

table(bwnet.4$colors)
table(bwnet.6$colors)
table(bwnet.8$colors)

# ============================================================
# 11. kME Summary Table
# ============================================================

library(WGCNA)

if (!exists("datExpr")) {
  vsd <- vst(dds)
  datExpr <- t(assay(vsd))
}

nets <- list(P4 = bwnet.4, P6 = bwnet.6, P8 = bwnet.8)

compute_kME_stats <- function(net, power_label, expr) {
  cat("\n>>> Computing kME for", power_label, "...\n")
  
  kme <- WGCNA::signedKME(expr, net$MEs)
  
  module_vec  <- net$colors
  self_colname <- paste0("kME", module_vec)
  self_kME <- kme[cbind(seq_len(nrow(kme)),
                        match(self_colname, colnames(kme)))]
  
  per_gene <- data.frame(
    gene     = rownames(kme),
    module   = module_vec,
    kME_self = as.numeric(self_kME),
    kme,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  
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

results <- list()
for (nm in names(nets)) {
  results[[nm]] <- compute_kME_stats(nets[[nm]], nm, datExpr)
}

kme_summary_all <- do.call(rbind, lapply(results, `[[`, "summary"))
rownames(kme_summary_all) <- NULL

write.csv(kme_summary_all,
          "mydata/myresults/kME_summary_all_powers.csv",
          row.names = FALSE)

for (nm in names(results)) {
  write.csv(results[[nm]]$per_gene,
            paste0("mydata/myresults/kME_pergene_", nm, ".csv"),
            row.names = FALSE)
}

for (nm in names(results)) {
  cat("\n===========================================\n")
  cat("Median kME —", nm, "\n")
  cat("===========================================\n")
  print(results[[nm]]$summary, row.names = FALSE)
}

modules_all <- sort(unique(kme_summary_all$module))
powers_all  <- sort(unique(kme_summary_all$power))

kme_wide <- data.frame(module = modules_all, stringsAsFactors = FALSE)
for (p in powers_all) {
  sub <- kme_summary_all[kme_summary_all$power == p, c("module", "median_kME")]
  kme_wide[[paste0("median_kME_", p)]] <-
    sub$median_kME[match(kme_wide$module, sub$module)]
}

kme_wide <- kme_wide[order(-kme_wide[[paste0("median_kME_", powers_all[length(powers_all)])]],
                           na.last = TRUE), ]
rownames(kme_wide) <- NULL

cat("\n===========================================\n")
cat("Median kME across powers\n")
cat("===========================================\n")
print(kme_wide, row.names = FALSE)

write.csv(kme_wide,
          "mydata/myresults/kME_median_wide.csv",
          row.names = FALSE)

# ---- gt table ----
library(gt)
library(scales)

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
      domain   = c(-0.5, 1),
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

# ============================================================
# 12. Hub-gene GO annotation tables (pink + cyan, power 8)
#     Produces Table A only: per-gene GO annotations for the
#     top 5 hub genes of each module, at P8.
# ============================================================

library(dplyr)
library(GO.db)
library(AnnotationDbi)

# ------------------------------------------------------------
# 12.1  Read and normalize the GO mapping
#        (skip if geneID2GO.raw already exists in your workspace)
# ------------------------------------------------------------
if (!exists("geneID2GO.raw")) {
  mappingFile <- "mydata/trinotate_annotation_GOblastx_forTopGO.txt"
  geneID2GO   <- readMappings(file = mappingFile)
  
  flatten <- function(x) {
    if (is.null(x))      return(character(0))
    if (is.character(x)) return(x)
    unique(unlist(x, use.names = FALSE, recursive = TRUE))
  }
  
  geneID2GO     <- lapply(geneID2GO, flatten)
  geneID2GO     <- geneID2GO[lengths(geneID2GO) > 0]
  geneID2GO.raw <- geneID2GO
  
  cat("Genes in GO mapping:", length(geneID2GO.raw), "\n")
}

# ------------------------------------------------------------
# 12.2  Shared helpers
# ------------------------------------------------------------
if (!exists("datExpr")) {
  vsd     <- vst(dds)
  datExpr <- t(assay(vsd))
}

safe_term <- function(g) {
  if (is.na(g) || !nzchar(g)) return(NA_character_)
  t <- tryCatch(AnnotationDbi::Term(g), error = function(e) NULL)
  if (is.null(t) || length(t) == 0 || is.na(t[1])) return(NA_character_)
  as.character(t[1])
}

safe_ont <- function(g) {
  if (is.na(g) || !nzchar(g)) return(NA_character_)
  o <- tryCatch(AnnotationDbi::Ontology(g), error = function(e) NULL)
  if (is.null(o) || length(o) == 0 || is.na(o[1])) return(NA_character_)
  as.character(o[1])
}

get_go <- function(id, mapping = geneID2GO.raw) {
  go <- mapping[[id]]
  if (is.null(go) || length(go) == 0) {
    return(data.frame(
      ID       = id,
      GO.ID    = NA_character_,
      Term     = NA_character_,
      Ontology = NA_character_,
      stringsAsFactors = FALSE
    ))
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

# ------------------------------------------------------------
# 12.3  Function: build the top-5 hub-gene GO table for one module
# ------------------------------------------------------------
hub_go_table <- function(net,
                         power_label,
                         module,
                         expr     = datExpr,
                         mapping  = geneID2GO.raw,
                         top_n    = 5) {
  
  # kME for every gene against every module
  kme       <- WGCNA::signedKME(expr, net$MEs)
  self_col  <- paste0("kME", module)
  
  if (!self_col %in% colnames(kme)) {
    stop("Module '", module, "' not found. Available: ",
         paste(sub("^kME", "", colnames(kme)), collapse = ", "))
  }
  
  hub.df <- data.frame(
    Gene   = colnames(expr),
    Module = net$colors,
    kME    = kme[[self_col]],
    stringsAsFactors = FALSE
  )
  
  # Filter to module FIRST, then sort
  hub.df <- subset(hub.df, Module == module)
  if (nrow(hub.df) == 0) {
    stop("No genes assigned to module '", module, "' in ", power_label, ".")
  }
  hub.df <- hub.df[order(-abs(hub.df$kME)), ]
  
  top    <- head(hub.df, min(top_n, nrow(hub.df)))
  top$ID <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", top$Gene)
  top$ID <- sub("^RINITY_", "TRINITY_", top$ID)
  
  # Long-format annotation table
  out        <- dplyr::bind_rows(lapply(top$ID, get_go, mapping = mapping))
  out$kME    <- top$kME[match(out$ID, top$ID)]
  out$Module <- module
  out$Power  <- power_label
  out        <- out[, c("Power", "Module", "ID", "GO.ID", "Term",
                        "Ontology", "kME")]
  
  # Write CSV (long format — Table A)
  outfile <- sprintf("mydata/myresults/top5_hub_GO_%s_%s.csv",
                     module, power_label)
  write.csv(out, outfile, row.names = FALSE)
  cat("Wrote:", outfile, "\n")
  
  as.data.frame(out)
}

# ------------------------------------------------------------
# 12.4  Pink module at P8
# ------------------------------------------------------------
top5_go_pink <- hub_go_table(bwnet.8, "P8", "pink", top_n = 5)

cat("\nTop 5 hub gene GO annotations — PINK, P8:\n")
print(top5_go_pink, row.names = FALSE)

# ------------------------------------------------------------
# 12.5  Cyan module at P8
# ------------------------------------------------------------
top5_go_cyan <- hub_go_table(bwnet.8, "P8", "cyan", top_n = 5)

cat("\nTop 5 hub gene GO annotations — CYAN, P8:\n")
print(top5_go_cyan, row.names = FALSE)
