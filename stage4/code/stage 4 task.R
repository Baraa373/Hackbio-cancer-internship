## load necessary libraries 

library("TCGAbiolinks")
library(SummarizedExperiment)
library(edgeR)
library(gplots)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(dplyr)
library(dplyr)
library(caret)  # for confusionMatrix
library(class)   # for KNN
library(Boruta)  # for feature selection
library(lattice)
# project information 
getProjectSummary("TCGA-LGG")

# download LGG dataset 
tcga_lgg <- GDCquery(project = "TCGA-LGG",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification")
GDCdownload(tcga_lgg) 
lgg_data <- GDCprepare(tcga_lgg) 
# explore IDH status in the metadata
lgg_data$barcode
lgg_data$paper_IDH.status
table(lgg_data$paper_IDH.status)
lgg_data$paper_IDH.specific.DNA.Methylation.Cluster # Lgm1-6
table(lgg_data$paper_IDH.specific.DNA.Methylation.Cluster)
lgg_data$paper_IDH.specific.RNA.Expression.Cluster # LGr1-4
table(lgg_data$paper_IDH.specific.RNA.Expression.Cluster)
# create a simple metadata 
meta <- data.frame("barcode" = lgg_data$barcode,
                   "IDH status" = lgg_data$paper_IDH.status)
View(meta)
table(is.na(meta)) # 21 NAs -> eliminate samples with no IDH status
meta <- na.omit(meta)
# get raw counts from lgg_data
lgg_rawdata <- assays(lgg_data) 
dim(lgg_rawdata$unstranded) 
View(lgg_rawdata$unstranded)

lgg <- lgg_rawdata$unstranded
samples <- meta$barcode
lgg <- lgg[, colnames(lgg) %in% samples]
all(colnames(lgg) %in% samples) # TRUE

#save the data set

# Transpose the gene expression data
lgg_t <- t(lgg)

# Convert to a data frame
lgg_df <- as.data.frame(lgg_t)

# Add metadata (like IDH status) as columns to the gene expression dataframe
lgg_merged <- merge(lgg_df, meta, by.x = "row.names", by.y = "barcode")

# View the merged data
View(lgg_merged)
str(lgg_merged)
# Normalize the data using log2 transformation
lgg_normalized <- log2(lgg_df + 1)

# Save the gene expression data to a CSV file
write.csv(lgg_df, file = "lgg_gene_expression.csv", row.names = TRUE)
# Save the metadata (barcodes and IDH statuses) to a CSV file
write.csv(meta, file = "lgg_metadata.csv", row.names = FALSE)

# Remove any non-numeric columns like 'IDH.status' if merged earlier
# Assuming row names are barcodes, we leave them as they are.
lgg_numeric <- lgg_df[, sapply(lgg_df, is.numeric)]

# Apply log2 transformation only to the numeric gene expression data
lgg_normalized <- log2(lgg_numeric + 1)

# Set the row names back to the barcodes
rownames(lgg_normalized) <- rownames(lgg_df)

# Check the first few rows to confirm
head(lgg_normalized)

# If you still need the 'IDH.status' or any metadata, merge it back afterward:
lgg_normalized_with_metadata <- cbind(lgg_normalized, lgg_merged$IDH.status)
colnames(lgg_normalized_with_metadata)[ncol(lgg_normalized_with_metadata)] <- "IDH.status"

# View the final data
View(lgg_normalized_with_metadata)

#Feature Selection
# Install Boruta if you haven't already
if (!requireNamespace("Boruta", quietly = TRUE)) {
  install.packages("Boruta")
}

library(Boruta)
# Ensure the IDH status is a factor
lgg_normalized_with_metadata$IDH.status <- as.factor(lgg_normalized_with_metadata$IDH.status)

# Prepare feature data (excluding the IDH status column)
feature_data <- lgg_normalized_with_metadata[, -ncol(lgg_normalized_with_metadata)]  # Exclude IDH status
idh_labels <- lgg_normalized_with_metadata$IDH.status  # IDH status

# Run Boruta
boruta_output <- Boruta(feature_data, idh_labels, doTrace = 2)

# Get selected features
selected_features <- getSelectedAttributes(boruta_output, withTentative = TRUE)

# Print selected features
print(selected_features)

# Split the data into training and testing sets (70% train, 30% test)
set.seed(234)  # For reproducibility
train_indices <- createDataPartition(lgg_normalized_with_metadata$IDH.status, p = 0.7, list = FALSE)

train_data <- lgg_normalized_with_metadata[train_indices, ]
test_data <- lgg_normalized_with_metadata[-train_indices, ]

# Set the features and labels for KNN
train_features <- train_data[, selected_features]
train_labels <- train_data$IDH.status
test_features <- test_data[, selected_features]
test_labels <- test_data$IDH.status

# Load the class library for KNN
library(class)

# Set the value of k
k_value <- 5  # You can adjust this value

# Run KNN
knn_predictions <- knn(train = train_features, test = test_features, cl = train_labels, k = k_value)

# Evaluate the model
library(caret)
confusion_matrix <- confusionMatrix(knn_predictions, test_labels)
print(confusion_matrix)

# Install randomForest if you haven't already
if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}

library(randomForest)

# Random Forest model training
rf_model <- randomForest(x = train_features, y = train_labels, importance = TRUE)

# Make predictions on the test set
rf_predictions <- predict(rf_model, newdata = test_features)

# Evaluate the Random Forest model
rf_confusion_matrix <- confusionMatrix(rf_predictions, test_labels)
print(rf_confusion_matrix)
# KNN Confusion Matrix Visualization with Custom Colors
ggplot(as.data.frame(confusion_matrix$table), aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "lightyellow", high = "darkorange") +  # Custom colors
  theme_minimal() +
  labs(title = "KNN Confusion Matrix", x = "Actual IDH Status", y = "Predicted IDH Status")













