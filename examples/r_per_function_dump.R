#!/usr/bin/env Rscript
# R driver for Notebook 3: dumps per-function outputs for R-Parity comparison.
suppressMessages({
  library(Statial)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript r_per_function_dump.R <fixture_path> <output_dir>")
}

fixture_path <- args[1]
output_dir   <- args[2]

input <- readRDS(fixture_path)
set.seed(42)

results <- list()

# getDistances
sce <- getDistances(input, maxDist = 200)
results$getDistances <- list(
  shape = dim(reducedDim(sce, "distances")),
  colnames = colnames(reducedDim(sce, "distances"))
)

# getAbundances
sce <- getAbundances(input, r = 200)
results$getAbundances <- list(
  shape = dim(reducedDim(sce, "abundances")),
  colnames = colnames(reducedDim(sce, "abundances"))
)

# Kontextual
kdf <- Kontextual(
  cells = input, r = 50,
  from = "Macrophages", to = "Keratin_Tumour",
  parent = c("Macrophages", "CD4_Cell"),
  image = "6", edgeCorrect = FALSE, window = "square"
)
results$Kontextual <- kdf

# parentCombinations
results$parentCombinations <- parentCombinations(
  all = c("tumour", "CD4", "CD8", "epithelial", "stromal"),
  tcells = c("CD4", "CD8"),
  tissue = c("epithelial", "stromal")
)

# makeWindow
cd <- as.data.frame(colData(input))
img6 <- cd[cd$imageID == "6", ]
results$makeWindow <- list(
  square = makeWindow(img6, window = "square"),
  convex = makeWindow(img6, window = "convex")
)

out_path <- file.path(output_dir, "r_per_function_outputs.json")
jsonlite::write_json(results, out_path, auto_unbox = TRUE, digits = NA, pretty = TRUE)
cat("Per-function outputs written to:", out_path, "\n")
