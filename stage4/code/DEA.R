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

# Install necessary packages if not installed
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) install.packages("EnhancedVolcano")

# Load libraries
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)

# Set the working directory (ensure correct escaping for backslashes)
setwd("D:\\bioinformatics\\hackbio")

# Import data with headers
upregulated_data <- read.csv("upregulated_genes_task5_updated.csv", header = TRUE)
downregulated_data <- read.csv("downregulated_genes_task5_updated.csv", header = TRUE)

# Remove rows with missing values (NA)
upregulated_data <- na.omit(upregulated_data)
downregulated_data <- na.omit(downregulated_data)

# Set the first column as row names and remove it from the data
rownames(upregulated_data) <- upregulated_data[, 1]
upregulated_data <- upregulated_data[, -1]  # Remove the first column from data

rownames(downregulated_data) <- downregulated_data[, 1]
downregulated_data <- downregulated_data[, -1]  # Remove the first column from data

# Combine upregulated and downregulated data
combined_data <- rbind(upregulated_data, downregulated_data)

# Calculate variance for each gene
gene_variances <- apply(combined_data, MARGIN = 1, FUN = var)

# Specify the number of top genes to keep (you can adjust this if necessary)
top_n <- 100

# Filter top genes based on variance
if (nrow(combined_data) < top_n) {
    top_genes <- combined_data  # Use the entire dataset if fewer genes than top_n
} else {
    top_genes_indices <- order(gene_variances, decreasing = TRUE)[1:top_n]
    top_genes <- combined_data[top_genes_indices, ]
}

# Check the dimensions to ensure you have the desired number of genes
dim(top_genes)

# Optional: Normalize the data (mean-center and scale)
normalized_data <- scale(top_genes)

# Create a heatmap of the top genes
pheatmap(normalized_data, 
         cluster_rows = TRUE,  # Cluster rows (genes)
         cluster_cols = TRUE,  # Cluster columns (samples)
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE,  # Show column names
         main = "Top 100 Genes by Variance")

# Optionally, you can save the heatmap to a file
# pheatmap(normalized_data, filename = "heatmap.png")

# Make a volcano plot using the combined dataset (if necessary)
EnhancedVolcano(
    combined_data,
    lab = rownames(combined_data),
    x = 'log2FoldChange',  # Replace with your actual column for log2 fold change
    y = 'padj',            # Replace with your actual column for adjusted p-value
    xlim = c(-6, 6),       # Adjust the x-axis limits as needed
    pCutoff = 0.05,        # P-value cutoff
    FCcutoff = 1,          # Fold-change cutoff
    pointSize = 3.0,
    labSize = 4.0,
    title = 'Volcano Plot',
    subtitle = 'Upregulated and Downregulated Genes',
    caption = 'log2FoldChange vs adjusted p-value'
)






























upregulated_counts <- normalized_counts[upregulated_genes, ]
downregulated_genes <- rownames(downregulated)
downregulated_counts <- normalized_counts[downregulated_genes, ]

