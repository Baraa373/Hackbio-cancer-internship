# Load necessary libraries
library(readr)  # For reading CSV files
library(ggplot2)  # For visualization

# Load the data
file_path <- "C:\\Users\\asus\\OneDrive\\Documents\\Table\\lgg_normalized_with_metadata.csv"
lgg_normalized_with_metadata <- read_csv(file_path)

# Check the structure of the dataset
str(lgg_normalized_with_metadata)
# Step 1: Extract only numeric columns for PCA
numeric_data <- lgg_normalized_with_metadata[, sapply(lgg_normalized_with_metadata, is.numeric)]

# Check the structure to ensure only numeric columns are included
str(numeric_data)
# Perform PCA on the numeric columns
# Extract only numeric columns from the dataset
numeric_data <- lgg_normalized_with_metadata %>% select(where(is.numeric))

# Perform PCA
pca_result <- prcomp(numeric_data, center = TRUE, scale. = TRUE)

# Identify columns with zero variance
zero_var_cols <- sapply(numeric_data, function(x) var(x) == 0)

# Get names of zero variance columns
constant_columns <- names(numeric_data)[zero_var_cols]
print(constant_columns)
# Remove constant columns
numeric_data_cleaned <- numeric_data[, !zero_var_cols]
# Perform PCA on the cleaned numeric columns
pca_result <- prcomp(numeric_data_cleaned, center = TRUE, scale. = TRUE)

# Create a data frame with PCA results and IDH status
pca_data <- data.frame(pca_result$x, IDH.status = lgg_normalized_with_metadata$IDH.status)

# Run K-means clustering on the first 10 PCA components
set.seed(123)  # For reproducibility
k_value <- 4  # Specify the number of clusters
kmeans_result <- kmeans(pca_data[, 1:10], centers = k_value, nstart = 25)

# Add cluster information to the PCA data
pca_data$cluster <- as.factor(kmeans_result$cluster)

# Cross-reference clusters with IDH status
cluster_idh_table <- table(pca_data$cluster, pca_data$IDH.status)
print(cluster_idh_table)

# Install ggplot2 if you haven't already
install.packages("ggplot2")

# Load ggplot2 for visualization
library(ggplot2)
# Choose one of the custom color palettes
custom_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#F0E442")

# Visualize the clusters with the chosen color palette
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster, shape = IDH.status)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering of Gliomas (PCA Reduced)", 
       x = "Principal Component 1", 
       y = "Principal Component 2") +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_minimal()
