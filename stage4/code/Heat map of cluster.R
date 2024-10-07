# Load necessary libraries
library(dplyr)

# Step 1: Load your gene expression data
gene_expression <- read.csv("C:\\Users\\asus\\OneDrive\\Documents\\lgg_gene_expression.csv", row.names = 1)

# Alternatively, select the top N genes by variance
top_n <- 100  # Specify the number of top genes you want to keep
if (nrow(filtered_genes) < top_n) {
  top_genes <- filtered_genes  # If filtered data has fewer genes than N, use it
} else {
  top_genes_indices <- order(gene_variances, decreasing = TRUE)[1:top_n]
  top_genes <- gene_expression[top_genes_indices, ]
}

# Check the dimensions of the top genes dataset
dim(top_genes)  # Should show the number of genes you want to keep

# Step 4: Subsetting Samples (if needed)
# Example: Randomly select a subset of samples
set.seed(42)  # For reproducibility
sample_size <- 50  # Specify the number of samples you want to keep
selected_samples <- sample(colnames(gene_expression), sample_size)

# Subset the gene expression data to keep only the selected samples
reduced_data <- top_genes[, selected_samples]

# Check the dimensions of the reduced dataset
dim(reduced_data)  # Should show reduced number of genes and samples

# Step 5: Save the reduced data if needed
write.csv(reduced_data, file = "C:\\Users\\asus\\OneDrive\\Documents\\lgg_reduced_gene_expression.csv", row.names = TRUE)
# Load necessary libraries
library(pheatmap)
library(dplyr)
library(tibble)
library(RColorBrewer)  # Load the RColorBrewer package for color palettes

# Step 1: Create clusters using k-means
num_clusters <- 2  # Set the number of clusters to 2
set.seed(42)
kmeans_result <- kmeans(t(reduced_data), centers = num_clusters)
cluster_assignments <- kmeans_result$cluster

# Step 2: Calculate cluster means
reduced_data_matrix <- as.matrix(reduced_data)
data_with_clusters <- data.frame(t(reduced_data_matrix), Cluster = factor(cluster_assignments))

# Calculate mean expression for each cluster
cluster_means <- data_with_clusters %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames(var = "Cluster")

# Step 3: Generate the heatmap with PuBu color palette
pheatmap(cluster_means, 
         main = "Heatmap of Cluster Means for Gene Expression Data",
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = brewer.pal(n = 9, name = "PuBu"),  # Using PuBu color palette
         fontsize_row = 10,       # Adjust row font size
         fontsize_col = 10,       # Adjust column font size
         fontsize = 12)           # Adjust title font size

# Save the heatmap as a PNG file (optional)
ggsave("C:\\Users\\asus\\OneDrive\\Documents\\heatmap_cluster_means.png")
