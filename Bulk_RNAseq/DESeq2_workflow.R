# DESeq2 Pipeline: Install packages once if needed
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

# DESeq2 Pipeline: Counts File Generation

library(DESeq2)

# Locate STAR-generated count files
file.list <- list.files(
  path = "./",
  pattern = "*ReadsPerGene.out.tab$"
)

# Load STAR count files, skipping header lines
counts.files <- lapply(file.list, read.table, skip = 4)

# Extract unstranded gene counts (column 2)
counts <- as.data.frame(
  sapply(counts.files, function(x) x[,2])
)

# Label samples using file names
colnames(counts) <- file.list

# Assign gene IDs as row names
row.names(counts) <- counts.files[[1]]$V1


# DESeq2 Pipeline: Create Metadata

# Define sample conditions
# Adjust order if file.list order changes
condition <- c(rep("AD",3), rep("Control",3))

# Create sample metadata table
sampleTable <- data.frame(
  sampleName = file.list,
  condition = condition
)


# Create DESeq2 dataset

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sampleTable,
  design = ~ condition
)


# DESeq2 Pipeline: Run Differential Expression Analysis

output <- DESeq(dds)

# Compare Alzheimer's disease vs Control
results_Control_AD <- results(
  output,
  contrast = c("condition","AD","Control")
)

# Sort results by adjusted p-value
results_Control_AD_PValue <- results_Control_AD[
  order(results_Control_AD$padj),
]

# View top results
head(results_Control_AD_PValue)

# Save results
write.csv(
  as.data.frame(results_Control_AD_PValue),
  file = "Alzheimers_DESeq2_results.csv"
)


# Convert ENSEMBL IDs to Gene Symbols

# Install once if needed:
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")

library(AnnotationDbi)
library(org.Hs.eg.db)

# Load DESeq2 results file
AD_Control <- read.table(
  "Control_AD_DE.csv",
  header = TRUE,
  sep = ','
)

# Extract ENSEMBL IDs
IDs <- c(AD_Control$X)

# Convert ENSEMBL IDs to gene symbols
AD_Control$Symbol <- mapIds(
  org.Hs.eg.db,
  IDs,
  'SYMBOL',
  'ENSEMBL'
)

# Set row names
rownames(AD_Control) <- AD_Control$X

# Save annotated results
write.csv(
  AD_Control,
  file = "Control_AD_DE_GeneIDs.csv"
)
