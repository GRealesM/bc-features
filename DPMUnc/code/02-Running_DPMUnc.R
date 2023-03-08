##################
# Running DPMUnc #
##################

# Date: 2022/08/24
# Guillermo Reales


# This script will run DPMUnc. I may adapt it to run in parallel. We'll see!

# Load packages
library(data.table)
#devtools::install_github("nichollskc/DPMUnc") # install latest version
library(DPMUnc)

# Helper functions
read_matrix <- function(x){
                  ds <- fread(x)
                  ds[, Trait:=NULL] # Remove Trait and label columns, keep just the numbers!
                  mat <- as.matrix(ds)
                  colnames(mat) <- NULL
                  return(mat)
}

# Load datasets
args <- commandArgs(trailingOnly = TRUE) # First argument should be the "experiment name", second the seed
exp <- args[1]
seed <- as.numeric(args[2])
delta <- paste0("../data/", exp, "_Delta.tsv")
var <- paste0("../data/", exp, "_Var.tsv")

message("Running DPMUnc on ", exp, " using seed = ", seed, ".")

obsData <- read_matrix(delta)
obsVars <- read_matrix(var)

# Standard run
DPMUnc(obsData, obsVars, saveFileDir = paste0("../results/", exp, "_", seed), seed = seed, kappa0 = 0.01, alpha0 = 2, beta0 = 0.1, nIts = 5000000, scaleData=TRUE)
message("Done!")

