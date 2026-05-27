#!/usr/bin/env Rscript
# R reference driver for py-statial parity tests.
# Loads the canonical fixture, runs Statial functions, writes JSON outputs.

suppressMessages({
  library(Statial)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript r_reference_driver.R <fixture_path> <output_path>")
}

fixture_path <- args[1]
output_path  <- args[2]

input <- readRDS(fixture_path)
set.seed(42)

results <- list()

# 1. getDistances
tryCatch({
  sce <- getDistances(input, maxDist = 200)
  dist <- as.matrix(reducedDim(sce, "distances"))
  # Replace NA with -1 for JSON compatibility
  dist[is.na(dist)] <- -1
  results$distances <- dist
  results$distances_colnames <- colnames(dist)
  results$distances_rownames <- rownames(dist)[1:100]  # first 100 for verification
}, error = function(e) {
  message("getDistances failed: ", e$message)
})

# 2. getAbundances
tryCatch({
  sce <- getAbundances(input, r = 200)
  abund <- as.matrix(reducedDim(sce, "abundances"))
  abund[is.na(abund)] <- -1
  results$abundances <- abund
  results$abundances_colnames <- colnames(abund)
}, error = function(e) {
  message("getAbundances failed: ", e$message)
})

# 3. Kontextual
tryCatch({
  kdf <- Kontextual(
    cells = input,
    r = 50,
    from = "Macrophages",
    to = "Keratin_Tumour",
    parent = c("Macrophages", "CD4_Cell"),
    image = "6",
    edgeCorrect = FALSE,
    window = "square"
  )
  results$kontextual <- kdf
}, error = function(e) {
  message("Kontextual failed: ", e$message)
})

# 4. calcContamination
tryCatch({
  sce <- calcContamination(input, num.trees = 50)
  contam <- as.data.frame(reducedDim(sce, "contaminations"))
  # Keep only numeric columns
  numeric_cols <- sapply(contam, is.numeric)
  contam_num <- contam[, numeric_cols, drop = FALSE]
  contam_mat <- as.matrix(contam_num)
  contam_mat[is.na(contam_mat)] <- -1
  results$contaminations <- contam_mat
  results$contaminations_colnames <- colnames(contam_mat)
}, error = function(e) {
  message("calcContamination failed: ", e$message)
})

# Serialize
jsonlite::write_json(results, output_path, auto_unbox = TRUE, digits = NA, pretty = TRUE)
cat("R reference outputs written to:", output_path, "\n")
