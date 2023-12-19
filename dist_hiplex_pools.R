library(writexl)
library(stringdist)
library(anticlust)
library(pheatmap)

# User Inputs (replace with your desired values)
input_file   <- "~/tmp/GMCL_pools.csv"
output_file  <- "~/tmp/GMCL_pools.pooled.csv"
K            <- 10               # number of pools
sgRNA_column <- "protospacer"    # Column name with the gRNA sequence
chr_column   <- "chromosome"     # Column name with the chromosome number
start_column <- "start"          # Column name with the start position

# Function to check chromosomal distances and exclude pairs within a specified distance
exclude_close_pairs <- function(df, max_distance) {
  num_gRNAs <- nrow(df)
  exclude_pairs <- matrix(FALSE, nrow=num_gRNAs, ncol=num_gRNAs)
  
  chr   <- df[, chr_column]
  start <- df[, start_column]

  # Iterate through gRNA pairs and check if gRNAs are on the same chromosome
  for (i in 1:(num_gRNAs - 1)) {
    for (j in (i + 1):num_gRNAs) {
      if (chr[i] == chr[j] & abs(start[i] - start[j]) < max_distance) {
        exclude_pairs[i, j] <- exclude_pairs[j, i] <- TRUE
      }
    }
  }
  
  return(exclude_pairs)
}

# Load the gRNA sequences from the CSV file
df <- read.csv(input_file)
distance_matrix <- as.matrix(stringdist::stringdistmatrix(df[, sgRNA_column], method="hamming"))

# Check chromosomal distances and exclude pairs within 10kb
distance_matrix[exclude_close_pairs(df, 1e4)] <- 0

# Distribute gRNAs using balanced distribution
set.seed(666)
df$pool <- anticlust::anticlustering(as.dist(distance_matrix), K, objective="distance")

# Print one distance matrix per pool
sapply(unique(df$pool), function(pool) {
  i   <- which(df$pool == pool)
  pal <- colorRampPalette(c("#F28E2B", "white", "#4E79A7"), bias=2)(max(c(10, nchar(df[, sgRNA_column]))))  # at least 10 breaks
  print(pheatmap::pheatmap(distance_matrix[i, i], color=pal, fontsize=6, main=paste("Pool", pool)))
})

# Output the result
l <- c(list(All_Pools=df), split(df, df$pool))
names(l) <- c("All_Pools", paste("Pool", names(l)[-1]))
writexl::write_xlsx(l, output_file)
