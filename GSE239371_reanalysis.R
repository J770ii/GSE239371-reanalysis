# ============================================================
# GSE239371 Reanalysis - Complete Pipeline
# Maternal microbiome programs developmental-stage-specific 
# DNA methylation in mouse hippocampus: a batch-aware reanalysis
#
# Author: Jalal Shegofte
# Date: 2026-10-15
# ============================================================

# ============================================================
# 0. SETUP
# ============================================================

# Set working directory
setwd("D:/GSE239371")

# Load libraries
library(limma)
library(org.Mm.eg.db)
library(GO.db)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(GEOquery)
library(R.utils)

# Create output folders
dir.create("objects", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# 1. DATA LOADING
# ============================================================

cat("=== Step 1: Loading GSE239371 data ===\n")

# Path to average beta file
beta_file <- "D:/GSE239371/GSE239371_average_beta.txt.gz/GSE239371_average_beta.txt.gz"

# Read the file
beta_df <- read.table(
    gzfile(beta_file),
    header = TRUE,
    sep = "\t",
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
)

cat("Dimensions:", dim(beta_df), "\n")

# Extract CpG IDs
cpg_ids <- beta_df$ID_REF

# ============================================================
# 2. QUALITY CONTROL AND FILTERING
# ============================================================

cat("\n=== Step 2: Quality control and filtering ===\n")

# Separate beta and detection p-values
col_names <- colnames(beta_df)
beta_cols <- grep("\\.AVG_Beta$", col_names, value = TRUE)
detp_cols <- grep("\\.Detection Pval$", col_names, value = TRUE)

beta_only <- beta_df[, beta_cols]
detp_only <- beta_df[, detp_cols]

# Convert to matrices
beta_matrix <- as.matrix(beta_only)
mode(beta_matrix) <- "numeric"
rownames(beta_matrix) <- make.unique(as.character(cpg_ids))

detP_matrix <- as.matrix(detp_only)
mode(detP_matrix) <- "numeric"
rownames(detP_matrix) <- make.unique(as.character(cpg_ids))

# Clean column names
colnames(beta_matrix) <- sub("\\.AVG_Beta$", "", colnames(beta_matrix))
colnames(detP_matrix) <- sub("\\.Detection Pval$", "", colnames(detP_matrix))

cat("Beta matrix:", dim(beta_matrix), "\n")

# Filter bad probes
detP_threshold <- 0.01
sample_fraction <- 0.05
bad_probes <- rowSums(detP_matrix > detP_threshold, na.rm = TRUE) > 
              (sample_fraction * ncol(detP_matrix))

beta_filtered <- beta_matrix[!bad_probes, ]

# Remove NA rows
na_rows <- rowSums(is.na(beta_filtered)) > 0
beta_clean <- beta_filtered[!na_rows, ]

cat("Clean beta matrix:", dim(beta_clean), "\n")

saveRDS(beta_clean, "objects/beta_clean.rds")
saveRDS(detP_matrix, "objects/detP_matrix.rds")

# ============================================================
# 3. METADATA
# ============================================================

cat("\n=== Step 3: Metadata preparation ===\n")

# Download metadata from GEO
options(timeout = 600)
gse <- getGEO("GSE239371", GSEMatrix = TRUE, getGPL = FALSE)
meta <- pData(gse[[1]])

meta_clean <- data.frame(
    GSM = rownames(meta),
    title = meta$title,
    group = meta[["experimental group:ch1"]],
    class = meta[["class:ch1"]],
    microbiota = meta[["gut microbiota:ch1"]],
    stringsAsFactors = FALSE
)

# Build final metadata
final_meta <- data.frame(
    col_index = 1:ncol(beta_clean),
    beta_name = colnames(beta_clean),
    GSM = meta_clean$GSM,
    title = meta_clean$title,
    group = meta_clean$group,
    class = meta_clean$class,
    microbiota = meta_clean$microbiota,
    stringsAsFactors = FALSE
)

saveRDS(final_meta, "objects/final_meta.rds")
write.csv(final_meta, "metadata/final_meta.csv", row.names = FALSE)

# ============================================================
# 4. BATCH ASSESSMENT
# ============================================================

cat("\n=== Step 4: Batch assessment ===\n")

# Extract Sentrix IDs from IDAT filenames
setwd("D:/GSE239371/IDATs_clean")
idat_files <- list.files(pattern = "_Grn\\.idat$")
sentrix <- sub(".*_(\\d+)_R.*", "\\1", idat_files)
position <- sub(".*_R(\\d+C\\d+)_.*", "\\1", idat_files)
gsm_id <- sub("_.*", "", idat_files)

batch_info <- data.frame(
    GSM = gsm_id,
    Sentrix = sentrix,
    Position = position,
    stringsAsFactors = FALSE
)

batch_info <- merge(batch_info, final_meta, by = "GSM")

# Cross-tabulation
batch_table <- table(batch_info$Sentrix, batch_info$group)
print(batch_table)

saveRDS(batch_info, "objects/batch_info.rds")
write.csv(as.data.frame.matrix(batch_table), 
          "results/tables/Table6_Batch_Effect.csv")

setwd("D:/GSE239371")

# ============================================================
# 5. DAM ANALYSIS (batch-clean)
# ============================================================

cat("\n=== Step 5: Dam differential methylation analysis ===\n")

dam_idx <- which(final_meta$class == "dam")
beta_dam <- beta_clean[, dam_idx]
meta_dam <- final_meta[dam_idx, ]

# M-value conversion
beta_dam <- beta_dam + 1e-6
M_dam <- log2(beta_dam / (1 - beta_dam))

# Design matrix
group_dam <- factor(meta_dam$group)
design_dam <- model.matrix(~ 0 + group_dam)
colnames(design_dam) <- levels(group_dam)

# Contrast
cont_dam <- makeContrasts(
    GM1_vs_GM4 = DamGM1 - DamGM4,
    levels = design_dam
)

# limma
fit_dam <- eBayes(contrasts.fit(lmFit(M_dam, design_dam), cont_dam), 
                   trend = TRUE)
res_dam <- topTable(fit_dam, coef = 1, number = Inf, adjust.method = "BH")

cat("Dam DMPs (FDR < 0.05):", sum(res_dam$adj.P.Val < 0.05), "\n")

saveRDS(res_dam, "results/dam_DMP_results.rds")
write.csv(res_dam, "results/tables/Table2_Dam_DMPs.csv")

# ============================================================
# 6. OFFSPRING ANALYSIS (batch-controlled)
# ============================================================

cat("\n=== Step 6: Offspring differential methylation analysis ===\n")

offspring_idx <- which(final_meta$class != "dam")
beta_off <- beta_clean[, offspring_idx]
meta_off <- final_meta[offspring_idx, ]

# Match batch info
batch_off <- batch_info$Sentrix[match(meta_off$GSM, batch_info$GSM)]

# M-value conversion
M_off <- log2((beta_off + 1e-6) / (1 - beta_off + 1e-6))

# Design with batch covariate
group <- factor(meta_off$group)
batch <- factor(batch_off)
design_off <- model.matrix(~ 0 + group + batch)
colnames(design_off) <- sub("^group", "", colnames(design_off))

# Contrasts
cont_off <- makeContrasts(
    Microbiota  = (GM1NonCF + GM1CF)/2 - (GM4NonCF + GM4CF)/2,
    CrossFoster = (GM1CF + GM4CF)/2 - (GM1NonCF + GM4NonCF)/2,
    Interaction = (GM1CF - GM1NonCF) - (GM4CF - GM4NonCF),
    levels = design_off
)

# limma
fit <- lmFit(M_off, design_off)
fit2 <- contrasts.fit(fit, cont_off)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# Extract results
results_off <- list()
for (ct in colnames(cont_off)) {
    res <- topTable(fit2, coef = ct, number = Inf, 
                    adjust.method = "BH", sort.by = "P")
    results_off[[ct]] <- res
    cat(ct, "DMPs (FDR<0.05):", sum(res$adj.P.Val < 0.05), "\n")
    write.csv(res, paste0("results/tables/DMP_", ct, ".csv"))
}

saveRDS(results_off, "objects/offspring_batch_controlled.rds")
saveRDS(fit2, "objects/limma_fit2_batch_controlled.rds")

# ============================================================
# 7. GENE ANNOTATION
# ============================================================

cat("\n=== Step 7: CpG annotation with TxDb ===\n")

# Load manifest
manifest_file <- "D:/GSE239371/IDATs_clean/GPL32685_MouseMethylation-12v1-0_A2.csv"
manifest <- read.csv(manifest_file, skip = 7, 
                     stringsAsFactors = FALSE, check.names = FALSE)

manifest_sub <- manifest[, c("IlmnID", "Name", "CHR", "MAPINFO")]
manifest_sub <- manifest_sub[!is.na(manifest_sub$CHR) & 
                              !is.na(manifest_sub$MAPINFO), ]

# Create GRanges
gr_cpg <- GRanges(
    seqnames = paste0("chr", manifest_sub$CHR),
    ranges = IRanges(start = as.numeric(manifest_sub$MAPINFO), width = 1),
    cpg_id = manifest_sub$Name
)

# Load TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- genes(txdb)

# Overlap
hits <- findOverlaps(gr_cpg, genes_gr)

cpg_gene <- data.frame(
    CpG = gr_cpg$cpg_id[queryHits(hits)],
    EntrezID = genes_gr$gene_id[subjectHits(hits)],
    CHR = as.character(seqnames(gr_cpg))[queryHits(hits)],
    Position = start(gr_cpg)[queryHits(hits)],
    stringsAsFactors = FALSE
)

cpg_gene$Gene_Symbol <- mapIds(org.Mm.eg.db,
                                keys = cpg_gene$EntrezID,
                                column = "SYMBOL",
                                keytype = "ENTREZID",
                                multiVals = "first")

saveRDS(cpg_gene, "objects/cpg_gene_map.rds")

# Annotate DMPs
# Dam
res_dam <- readRDS("results/dam_DMP_results.rds")
sig_dam <- res_dam[res_dam$adj.P.Val < 0.05, ]
sig_dam$CpG <- rownames(sig_dam)
sig_dam_annot <- merge(sig_dam, cpg_gene, by = "CpG", all.x = TRUE)
genes_dam <- unique(na.omit(sig_dam_annot$Gene_Symbol))
saveRDS(genes_dam, "objects/dam_genes.rds")
write.csv(sig_dam_annot, "results/dam_DMP_annotated.csv", row.names = FALSE)

# Offspring
res_off <- results_off$Microbiota
sig_off <- res_off[res_off$adj.P.Val < 0.05, ]
sig_off$CpG <- rownames(sig_off)
sig_off_annot <- merge(sig_off, cpg_gene, by = "CpG", all.x = TRUE)
genes_off <- unique(na.omit(sig_off_annot$Gene_Symbol))
saveRDS(genes_off, "objects/off_genes.rds")
write.csv(sig_off_annot, "results/off_DMP_annotated.csv", row.names = FALSE)

cat("Dam genes:", length(genes_dam), "\n")
cat("Offspring genes:", length(genes_off), "\n")
cat("Common genes:", length(intersect(genes_dam, genes_off)), "\n")

# ============================================================
# 8. GO ENRICHMENT (custom hypergeometric test)
# ============================================================

cat("\n=== Step 8: GO enrichment ===\n")

run_go_manual <- function(query_genes, ontology = "BP", 
                           p_cutoff = 0.05, min_genes = 3) {
    
    query_entrez_df <- suppressMessages(
        select(org.Mm.eg.db, keys = query_genes, 
               columns = "ENTREZID", keytype = "SYMBOL")
    )
    query_entrez_df <- query_entrez_df[!is.na(query_entrez_df$ENTREZID), ]
    query_ids <- unique(query_entrez_df$ENTREZID)
    
    all_entrez <- keys(org.Mm.eg.db, keytype = "ENTREZID")
    bg_go <- suppressMessages(
        select(org.Mm.eg.db, keys = all_entrez,
               columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID")
    )
    bg_go <- bg_go[!is.na(bg_go$GO) & bg_go$ONTOLOGY == ontology, ]
    
    N <- length(unique(bg_go$ENTREZID))
    K <- length(query_ids)
    bg_by_go <- split(bg_go$ENTREZID, bg_go$GO)
    
    results_list <- lapply(names(bg_by_go), function(go_term) {
        genes_in_bg <- unique(bg_by_go[[go_term]])
        m <- length(genes_in_bg)
        if (m < min_genes) return(NULL)
        genes_in_q <- intersect(query_ids, genes_in_bg)
        q <- length(genes_in_q)
        if (q == 0) return(NULL)
        p <- phyper(q - 1, m, N - m, K, lower.tail = FALSE)
        data.frame(GO_ID = go_term, Count = q, SetSize = m, pvalue = p,
                   geneID_entrez = paste(genes_in_q, collapse = "/"),
                   stringsAsFactors = FALSE)
    })
    
    results <- do.call(rbind, results_list)
    if (is.null(results) || nrow(results) == 0) return(NULL)
    
    results$p.adjust <- p.adjust(results$pvalue, method = "BH")
    go_desc <- suppressMessages(
        select(GO.db, keys = results$GO_ID, columns = "TERM", keytype = "GOID")
    )
    results$Description <- go_desc$TERM[match(results$GO_ID, go_desc$GOID)]
    results <- results[order(results$pvalue), ]
    results <- results[results$p.adjust < p_cutoff, ]
    
    if (nrow(results) > 0) {
        results$geneID <- sapply(results$geneID_entrez, function(x) {
            ids <- strsplit(x, "/")[[1]]
            syms <- suppressMessages(
                select(org.Mm.eg.db, keys = ids, 
                       columns = "SYMBOL", keytype = "ENTREZID")
            )
            paste(na.omit(syms$SYMBOL), collapse = "/")
        })
        results$geneID_entrez <- NULL
    }
    return(results)
}

# Run for dams
go_dam <- run_go_manual(genes_dam, ontology = "BP", p_cutoff = 0.05)
if (!is.null(go_dam)) {
    write.csv(go_dam, "results/tables/Table5_GO_Dams.csv", row.names = FALSE)
    saveRDS(go_dam, "objects/go_dam_manual.rds")
}

# Run for offspring
go_off <- run_go_manual(genes_off, ontology = "BP", p_cutoff = 0.05)
if (!is.null(go_off)) {
    write.csv(go_off, "results/tables/Table4_GO_Offspring.csv", row.names = FALSE)
    saveRDS(go_off, "objects/go_off_manual.rds")
}

# ============================================================
# 9. FIGURES
# ============================================================

cat("\n=== Step 9: Figures ===\n")

# Font setup
windowsFonts(Times = windowsFont("Times New Roman"))

theme_times <- function(base_size = 11) {
    theme_bw(base_size = base_size, base_family = "Times") +
    theme(
        text = element_text(family = "Times", size = 11, face = "bold"),
        plot.title = element_text(family = "Times", size = 11, face = "bold", hjust = 0.5),
        axis.title = element_text(family = "Times", size = 11, face = "bold"),
        axis.text = element_text(family = "Times", size = 11, face = "bold"),
        legend.title = element_text(family = "Times", size = 11, face = "bold"),
        legend.text = element_text(family = "Times", size = 11, face = "bold"),
        panel.grid = element_blank()
    )
}

# --- Figure 1: PCA ---
pca <- prcomp(t(M_off), scale. = TRUE)
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    Group = factor(meta_off$group,
                   levels = c("GM1CF", "GM1NonCF", "GM4CF", "GM4NonCF"))
)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(level = 0.68, linetype = 2, alpha = 0.4) +
    scale_color_manual(values = c(
        "GM1CF" = "#E41A1C", "GM1NonCF" = "#377EB8",
        "GM4CF" = "#4DAF4A", "GM4NonCF" = "#984EA3")) +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    labs(title = "PCA of offspring hippocampal DNA methylation",
         x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)")) +
    theme_times()

ggsave("figures/Figure1_PCA.png", p1, width = 8, height = 6, dpi = 300)

# --- Figure 2: Volcano ---
res <- results_off$Microbiota
res <- res[!is.na(res$logFC) & !is.na(res$adj.P.Val), ]
res$logFC_plot <- pmax(pmin(res$logFC, 3), -3)

volcano_df <- data.frame(
    CpG = rownames(res),
    logFC = res$logFC_plot,
    neglogP = -log10(res$adj.P.Val)
)
volcano_df$Sig <- "Not significant"
volcano_df$Sig[volcano_df$logFC > 0.05 & volcano_df$neglogP > -log10(0.05)] <- "Hypermethylated"
volcano_df$Sig[volcano_df$logFC < -0.05 & volcano_df$neglogP > -log10(0.05)] <- "Hypomethylated"
volcano_df$Sig <- factor(volcano_df$Sig,
                         levels = c("Hypermethylated", "Hypomethylated", "Not significant"))

top_cpg <- volcano_df[volcano_df$Sig != "Not significant", ]
top_cpg <- top_cpg[order(top_cpg$neglogP, decreasing = TRUE), ][1:min(10, nrow(top_cpg)), ]

p2 <- ggplot(volcano_df, aes(x = logFC, y = neglogP, color = Sig)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(
        "Hypermethylated" = "#E41A1C",
        "Hypomethylated" = "#377EB8",
        "Not significant" = "grey75")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", linewidth = 0.4) +
    geom_text_repel(data = top_cpg, aes(label = CpG),
                    size = 3, family = "Times", fontface = "bold") +
    labs(title = "Differentially methylated positions (Microbiota)",
         x = expression(Delta * beta ~ "(M-value difference, clipped to ±3)"),
         y = expression(-log[10] ~ "(adjusted P-value)"),
         color = "Significance") +
    theme_times()

ggsave("figures/Figure2_Volcano.png", p2, width = 8, height = 6, dpi = 300)

# --- Figure 3: Heatmap ---
top50 <- head(rownames(res)[order(res$adj.P.Val)], 50)
mat_top <- beta_off[top50, ]

annotation_col <- data.frame(
    Group = factor(meta_off$group,
                   levels = c("GM1CF", "GM1NonCF", "GM4CF", "GM4NonCF"))
)
rownames(annotation_col) <- colnames(mat_top)

png("figures/Figure3_Heatmap.png", width = 10, height = 10, units = "in", res = 300)
pheatmap(mat_top, annotation_col = annotation_col,
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 11, main = "Top 50 Differentially Methylated CpGs",
         border_color = NA)
dev.off()

# --- Figure 4: GO Barplot ---
if (!is.null(go_off) && nrow(go_off) > 0) {
    go_top <- head(go_off, 15)
    go_top$Description <- factor(go_top$Description,
                                  levels = rev(go_top$Description))
    
    p4 <- ggplot(go_top, aes(x = Count, y = Description, 
                              fill = -log10(p.adjust))) +
        geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
        scale_fill_gradient(low = "#FEE0B6", high = "#E41A1C",
                            name = expression(-log[10] ~ "(FDR)")) +
        labs(title = "Gene Ontology enrichment — Offspring DMPs",
             x = "Gene count", y = NULL) +
        theme_times()
    
    ggsave("figures/Figure4_GO_barplot.png", p4, width = 9, height = 7, dpi = 300)
}

# ============================================================
# 10. SESSION INFO
# ============================================================

cat("\n=== Analysis complete ===\n")
sink("results/session_info.txt")
cat("Session Info - GSE239371 Reanalysis\n")
cat("Date:", Sys.time(), "\n\n")
sessionInfo()
sink()

cat("All results saved in: D:/GSE239371/\n")