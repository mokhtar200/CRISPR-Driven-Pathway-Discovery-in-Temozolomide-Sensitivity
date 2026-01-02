#======================#
# 1. Load Dependencies #
#======================#
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ReactomePA"))
install.packages(c("ggplot2", "dplyr", "readxl"))

# Load Required Libraries
library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
#======================#
# 2. Load & Prepare Data
#======================#
# 1. Reconfirm column names are correct in screen_data
screen_data <- read_excel("D:/CRISPR/mmc4.xlsx", sheet = 1)
screen_data <- as_tibble(screen_data)
colnames(screen_data)[1] <- "gene"

# 2. Clean and prepare the full gene list
full_gene_list <- screen_data %>%
  dplyr::select(gene, G510_deltaZ) %>%
  dplyr::filter(!is.na(G510_deltaZ)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# 3. Get top hits (Z ≥ 2.0)
hits <- full_gene_list %>%
  dplyr::filter(G510_deltaZ >= 2.0) %>%
  dplyr::arrange(desc(G510_deltaZ)) %>%
  dplyr::rename(score = G510_deltaZ)

#==========================#
# 4. Visualization: Top Hits
#==========================#
hits %>%
  top_n(20, score) %>%
  ggplot(aes(x = reorder(gene, score), y = score)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  labs(title = "Top 20 CRISPR Hits (Temozolomide Sensitivity)",
       x = "Gene", y = "Delta Z-score") +
  theme_minimal()

#==========================#
# Functional Enrichment Analysis (GO & KEGG)
#==========================#
# Gene Set Enrichment Analysis - GO BP
gsea_go <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pAdjustMethod= "BH",
  keyType      = "ENTREZID",
  verbose      = FALSE
)

# Plot top enriched pathways
if (nrow(as.data.frame(gsea_go)) > 0) {
  ridgeplot(gsea_go, showCategory = 10, fill = "p.adjust") +
    ggtitle("GSEA - GO Biological Processes")
} else {
  message("No significant GSEA-GO results.")
}

write.csv(as.data.frame(gsea_go), "GSEA_GO_BP_results.csv", row.names = FALSE)


#======================#
# 8. Highlight Key Genes (e.g. MGMT, BARD1)
#======================#
highlight_genes <- c("MGMT", "BARD1", "BCL6", "SORBS2", "FAM179A", "MRPS35")
highlight_ids <- bitr(highlight_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Find pathways where they occur
reactome_hits <- gsea_reactome@result %>%
  filter(grepl(paste(highlight_ids$ENTREZID, collapse = "|"), core_enrichment))

print(reactome_hits[, c("ID", "Description", "core_enrichment")])


# =====================================================
#Volcano Plot for Whole Screen: CRISPR Knockout Effect
#======================================================

# Step 1: Prepare volcano data
volcano_data <- screen_data %>%
  mutate(
    significant = case_when(
      G510_deltaZ >= 2  ~ "Up",
      G510_deltaZ <= -2 ~ "Down",
      TRUE              ~ "NS"
    ),
    logZ = -log10(abs(G510_deltaZ))
  )

# Step 2: Select top 10 up and top 10 down genes
label_genes <- volcano_data %>%
  filter(significant != "NS") %>%
  group_by(significant) %>%
  slice_max(order_by = abs(G510_deltaZ), n = 10)

# Step 3: Volcano plot with selective labels
ggplot(volcano_data, aes(x = G510_deltaZ, y = logZ, color = significant)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3.5,
    max.overlaps = 100,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("firebrick", "grey", "dodgerblue")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Top 10 Up & Down-Regulated Genes",
    x = "Delta Z Score",
    y = "-log10(abs(Z))"
  )


# ===================================================
# MGMT Methylation-Specific Comparison: Split data by MGMT methylation status and compare which genes are selectively sensitizing in MGMT-methylated cells (e.g., Av_MGMT_methylated column).
#======================================

# Subset methylated-specific scores
methylated <- screen_data %>%
  dplyr::select(gene, Av_MGMT_methylated) %>%
  filter(Av_MGMT_methylated >= 2)

# Plot
ggplot(methylated, aes(x = reorder(gene, Av_MGMT_methylated), y = Av_MGMT_methylated)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "MGMT Methylated: Top TMZ Sensitizers",
       x = "Gene", y = "Z-score") +
  theme_minimal()

#========================================
# Cell Line-Specific Sensitivity Profiling
# Compare which genes are commonly or uniquely sensitizing across multiple glioblastoma lines (G361, G510, G523...).
#================================================================================================

# Get sets
g361_hits <- screen_data %>% filter(G361_deltaZ >= 2) %>% pull(gene)
g510_hits <- screen_data %>% filter(G510_deltaZ >= 2) %>% pull(gene)
g523_hits <- screen_data %>% filter(G523_deltaZ >= 2) %>% pull(gene)

# Venn or upset plot
library(VennDiagram)
venn.plot <- venn.diagram(list(G361 = g361_hits, G510 = g510_hits, G523 = g523_hits),
                          filename = NULL)
grid.draw(venn.plot)

#================================
# Gene Set Variation Analysis (GSVA): Assess pathway activity patterns across all cell lines in your screen.
#================================
library(GSVA)

expr_matrix <- as.matrix(screen_data[, c("G361_deltaZ", "G510_deltaZ", "G523_deltaZ")])
rownames(expr_matrix) <- screen_data$gene

msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)

gsva_res <- gsva(expr_matrix, hallmark_list, method = "gsva")
heatmap(gsva_res)


#==============================
#
#===================================

library(depmap)
 
# Get gene effect data
crispr_data <- depmap_crispr()  # This downloads the full CRISPR gene effect dataset

# Filter for MGMT scores across all cell lines
mgmt_scores <- crispr_data %>%
  filter(gene_name == "MGMT") %>%
  dplyr::select(depmap_id, mgmt_score = dependency)

# Merge with full data
merged <- crispr_data %>%
  filter(gene_name != "MGMT") %>%
  inner_join(mgmt_scores, by = "depmap_id")

# Correlation of all genes with MGMT
co_dep <- merged %>%
  group_by(gene_name) %>%
  summarise(correlation = cor(dependency, mgmt_score, method = "pearson"),
            n = n()) %>%
  arrange(desc(abs(correlation)))

# View top co-dependencies
head(co_dep, 10)


library(ggplot2)

ggplot(co_dep %>% top_n(20, abs(correlation)),
       aes(x = reorder(gene_name, correlation), y = correlation, fill = correlation > 0)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top 20 Gene Co-Dependencies with MGMT",
       x = "Gene", y = "Correlation (Pearson)") +
  theme_minimal()

#=================================

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

brain_genes <- getBM(attributes = c("hgnc_symbol", "description"),
                     filters = "hgnc_symbol",
                     values = hits$gene,
                     mart = ensembl)

# Example: suppose brain_genes$hgnc_symbol contains known brain-expressed genes
hits$brain_hit <- ifelse(hits$gene %in% brain_genes$hgnc_symbol, "Brain-Expressed", "Other")

library(ggplot2)
ggplot(hits %>% top_n(20, score), aes(x = reorder(gene, score), y = score, fill = brain_hit)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Brain-Expressed" = "blue", "Other" = "gray")) +
  labs(title = "Top CRISPR Hits with Brain-Expressed Genes Highlighted",
       x = "Gene", y = "Delta Z-Score") +
  theme_minimal()

#==============================
#  Correlation Between MGMT Methylation and CRISPR Sensitivity
#==========================
library(ggpubr)

screen_data %>%
  ggplot(aes(x = Av_MGMT_methylated, y = G510_deltaZ)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  stat_cor(method = "pearson", label.x = 1, label.y = 3) +
  labs(title = "Correlation: MGMT Methylation vs. CRISPR Sensitivity",
       x = "Average MGMT Methylation",
       y = "CRISPR Z-score (G510)") +
  theme_minimal()


#=========================================
#Heatmap of Top 50 Hits Across Cell Lines
#=========================================
library(tibble)
library(pheatmap)

top_hits <- screen_data %>%
  arrange(desc(G510_deltaZ)) %>%
  head(50) %>%
  column_to_rownames("gene")

expr_mat <- as.matrix(top_hits[, c("G361_deltaZ", "G510_deltaZ", "G523_deltaZ")])

pheatmap(expr_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Top 50 CRISPR Hits Across Cell Lines",
         color = colorRampPalette(c("navy", "white", "firebrick"))(50))

# ============================
# Dependency Score Comparison Between Methylated vs. Unmethylated
# ============================
screen_data %>%
  mutate(MGMT_group = ifelse(Av_MGMT_methylated >= 2, "Methylated", "Unmethylated")) %>%
  ggplot(aes(x = MGMT_group, y = G510_deltaZ, fill = MGMT_group)) +
  geom_boxplot() +
  labs(title = "CRISPR Sensitivity by MGMT Methylation",
       y = "Delta Z-score (G510)", x = "MGMT Methylation Status") +
  theme_minimal() +
  scale_fill_manual(values = c("dodgerblue", "firebrick"))



#==========================================
# DNA Repair Genes In CRISPR Knockout
#==========================================
# Define DNA repair genes
repair_genes <- c("BRCA1", "BRCA2", "RAD51", "XRCC1", "MLH1", "MSH2")

# Create annotated data
annotated_data <- screen_data %>%
  mutate(significant = ifelse(G510_deltaZ >= 2, "Up",
                              ifelse(G510_deltaZ <= -2, "Down", "NS")),
         is_repair = ifelse(gene %in% repair_genes, TRUE, FALSE),
         logZ = -log10(abs(G510_deltaZ)))

# Plot with labels for DNA repair genes
ggplot(annotated_data, aes(x = G510_deltaZ, y = logZ, color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = subset(annotated_data, is_repair == TRUE),
                  aes(label = gene), size = 3.5, color = "black") +
  labs(title = "DNA Repair Genes in CRISPR Knockout Volcano",
       x = "Delta Z Score", y = "-log10(abs(Z))") +
  scale_color_manual(values = c("firebrick", "grey", "dodgerblue")) +
  theme_minimal()

#================================================
# Top 10 Sensitive and Resistant Genes
#================================================
# Select Top 10 Sensitive (highest positive delta Z)
top_sensitive <- screen_data %>%
  filter(G510_deltaZ >= 2) %>%
  arrange(desc(G510_deltaZ)) %>%
  head(10)

# Select Top 10 Resistant (most negative delta Z)
top_resistant <- screen_data %>%
  filter(G510_deltaZ <= -2) %>%
  arrange(G510_deltaZ) %>%
  head(10)

# Combine
top_combined <- rbind(top_sensitive, top_resistant)

# Create expression matrix
clustering_mat <- as.matrix(top_combined[, c("G361_deltaZ", "G510_deltaZ", "G523_deltaZ")])
rownames(clustering_mat) <- top_combined$gene

# Create annotation row
annotation_df <- data.frame(
  Type = ifelse(top_combined$G510_deltaZ >= 2, "Sensitive", "Resistant")
)
rownames(annotation_df) <- top_combined$gene

# Heatmap
pheatmap(clustering_mat,
         annotation_row = annotation_df,
         main = "Top 10 Sensitive and Resistant Genes",
         fontsize_row = 10,
         cluster_cols = TRUE,
         cluster_rows = TRUE)


#======================#
# 9. Export Results
#======================#
write.csv(hits, "CRISPR_TEMO_TopHits.csv", row.names = FALSE)
write.csv(as.data.frame(ego_bp), "GO_BP_ORA.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "KEGG_ORA.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_go), "GO_BP_GSEA.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), "KEGG_GSEA.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_reactome), "Reactome_GSEA.csv", row.names = FALSE)
