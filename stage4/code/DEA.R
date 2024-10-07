# Install necessary packages
install.packages(c("ggplot2", "dplyr", "pheatmap", "RColorBrewer", "ggrepel"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Load libraries
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("D:\\bioinformatics\\hackbio")

# Load count data and sample information
count_data <- read.csv('Igg_gene_expression_updated.csv', header = TRUE, row.names = 1)
sample_info <- read.csv('lgg_metadata.csv', header = TRUE, stringsAsFactors = TRUE, fill = TRUE)
# Check the structure of the sample_info to confirm column names
head(sample_info)

# If needed, adjust the column names here
# Ensure that SampleID column exists in sample_info
sample_info$Sample.Type<- gsub("-", ".", sample_info$Sample.Type)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~Sample.Type)

# Ensure Sample.Type is a factor with specified levels
dds$Sample.Type <- factor(dds$Sample.Type, levels = c("WT", "Mutant"))

# Perform DESeq analysis
dds <- DESeq(dds)

# Extract results
deseq_results <- results(dds)

# Convert results to a data frame
deseq_results <- as.data.frame(deseq_results)

# Filter significant results
filtered <- deseq_results %>% filter(padj < 0.05)
filtered <- filtered %>% filter(abs(log2FoldChange) > 1)

# Identify upregulated and downregulated genes
upregulated <- deseq_results[deseq_results$log2FoldChange > 1 & deseq_results$padj < 0.05, ]
downregulated <- deseq_results[deseq_results$log2FoldChange < -1 & deseq_results$padj < 0.05, ]

upregu# Save upregulated and downregulated genes to CSV files
write.csv(upregulated, file = "upregulated_genes.csv", row.names = TRUE)
write.csv(downregulated, file = "downregulated_genes.csv", row.names = TRUE)

# Extract normalized count data
normalized_counts <- counts(dds, normalized = TRUE)

# Subset normalized counts for upregulated and downregulated genes
upregulated_genes <- rownames(upregulated)
upregulated_counts <- normalized_counts[upregulated_genes, ]
downregulated_genes <- rownames(downregulated)
downregulated_counts <- normalized_counts[downregulated_genes, ]

# Create heatmap for upregulated genes
pheatmap(upregulated_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         main = "Heatmap of Upregulated Genes")

# Create heatmap for downregulated genes
pheatmap(downregulated_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         main = "Heatmap of Downregulated Genes")

# Apply variance stabilizing transformation (VST)
vst_data <- vst(dds, blind = FALSE)

# Save the upregulated genes heatmap
png("upregulated_heatmap.png", width = 800, height = 800)
pheatmap(assay(vst_data)[rownames(upregulated_genes), ], cluster_rows = TRUE, show_rownames = TRUE, 
         cluster_cols = TRUE, annotation_col = as.data.frame(colData(dds)[, "Sample.Type"]))
dev.off()

# Load EnhancedVolcano library
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# Create Volcano plot using DESeq results
EnhancedVolcano(deseq_results,
                lab = rownames(deseq_results),
                x = 'log2FoldChange',  # x-axis for fold change
                y = 'pvalue',          # y-axis for p-values
                pCutoff = 0.05,        # Adjust as needed
                FCcutoff = 1,          # Fold-change threshold
                title = 'Volcano Plot',
                xlab = 'Log2 Fold Change',
                ylab = '-Log10 p-value',
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendLabels = c('NS', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                labSize = 3.0)

# Add a column to categorize genes as significant or not
deseq_results$significance <- ifelse(deseq_results$pvalue < 0.05 & abs(deseq_results$log2FoldChange) >= 1, 
                                     "Significant", "Not Significant")

# Create a volcano plot using ggplot2
ggplot(deseq_results, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("red", "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  theme(legend.position = "right")
