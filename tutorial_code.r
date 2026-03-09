rm(list = ls())
source("comp_function.R")
source("seminr.R")

library(readxl)
library(dplyr)
library(psych)
library(caret)
library(Metrics)
library(lmtest)
library(MASS)

# ==============================================================================
# Step 0: Data Preprocessing
#
# Read the raw dataset from Excel. The first 30 columns contain the raw
# compositional measurements for Attitude (cols 1-12), PBC (cols 13-21),
# and BI (cols 22-30). Each construct has 3-part compositions repeated
# multiple times (e.g., Attitude has 4 compositions × 3 parts = 12 cols).
# Columns 31-32 are the already ILR-transformed Subject Norms (SNnew).
# Zero values are replaced with 0.001 to avoid log(0) in ILR transformation.
# ==============================================================================

data <- read_excel("data_sub.xlsx", sheet = 1)
SN   <- data[, 31:32]
data <- data[, 1:30]

colnames(data) <- c(
  multi_items("Attitude", 1:12),
  multi_items("PBC", 1:9),
  multi_items("BI", 1:9)
)

data <- as.matrix(data)
data[data == 0] <- 0.001

# ==============================================================================
# Step 1: ILR Transformation of Compositional Variables
#
# Each group of 3 columns represents one 3-part composition. There are 10
# such groups across the 30 columns (Attitude: 4 groups, PBC: 3, BI: 3).
#
# 1) Apply compositional closure so each row sums to 1 within its 3-part group.
# 2) Apply the ILR (Isometric Log-Ratio) transformation to map each 3-part
#    composition to 2 real-valued coordinates in the Aitchison simplex.
#    This reduces 30 columns to 20, plus the 2 SNnew columns = 22 total.
# ==============================================================================

becomp <- function(i) Closure(as.data.frame(data)[, (3 * i - 2):(3 * i)])
do.ilr <- function(i) ilr(data[, (3 * i - 2):(3 * i)])

data     <- do.call("cbind", apply(matrix(1:10), 1, becomp))
data.ilr <- do.call("cbind", lapply(matrix(1:10), do.ilr))

# After ILR: Attitude has 8 cols (4 compositions × 2), PBC 6, BI 6
colnames(data.ilr) <- c(
  multi_items("Attitude", 1:8),
  multi_items("PBC", 1:6),
  multi_items("BI", 1:6)
)

# Combine ILR-transformed indicators with the pre-transformed SNnew columns
data.ilr <- cbind(data.ilr, SN)
obsData  <- data.ilr

head(obsData)
dim(obsData)

# ==============================================================================
# Step 2: Specify Measurement and Structural Models
#
# Measurement model (outer model): defines which ILR-transformed indicators
# belong to which construct. All constructs use Mode A (correlation weights).
#   - Attitude: 8 indicators (Attitude1–8)
#   - BI:       6 indicators (BI1–6)
#   - PBC:      6 indicators (PBC1–6)
#   - SNnew:    2 indicators (SNnew1–2)
#
# Structural model (inner model): defines the path relationships between
# constructs following the Theory of Planned Behavior (TPB):
#   SNnew → Attitude, SNnew → PBC, SNnew → BI, Attitude → BI, PBC → BI
# ==============================================================================

mm_model <- constructs(
  composite("Attitude", multi_items("Attitude", 1:8), weights = mode_A),
  composite("BI",       multi_items("BI", 1:6),       weights = mode_A),
  composite("PBC",      multi_items("PBC", 1:6),      weights = mode_A),
  composite("SNnew",    multi_items("SNnew", 1:2),     weights = mode_A)
)

sm_model <- relationships(
  paths(from = "SNnew",    to = "Attitude"),
  paths(from = "SNnew",    to = "PBC"),
  paths(from = "SNnew",    to = "BI"),
  paths(from = "Attitude", to = "BI"),
  paths(from = "PBC",      to = "BI")
)

# Convert model specifications into matrix format required by the PLS algorithm
smMatrix <- sm_model
measurement_model <- mm_model

recognized_constructs <- c("composite", "reflective", "higher_order_composite")
construct_measurements <- measurement_model[names(measurement_model) %in% recognized_constructs]
mmMatrix <- matrix(
  unlist(construct_measurements),
  ncol = 3, byrow = TRUE,
  dimnames = list(NULL, c("construct", "measurement", "type"))
)

constructs <- construct_names(smMatrix)
mmVariables <- mmMatrix[mmMatrix[, "construct"] %in% constructs, "measurement"]

# Determine the weighting mode (A or B) for each construct
measurement_mode_scheme <- sapply(
  unique(c(sm_model[, 1], sm_model[, 2])),
  get_measure_mode, mmMatrix, USE.NAMES = TRUE
)

# ==============================================================================
# Step 3: Standardize the Data
#
# Center and scale all measurement variables to zero mean and unit variance.
# This is required before running the PLS-SEM algorithm. The mean and SD are
# saved for later denormalization in the prediction step.
# ==============================================================================

normData <- standardize_safely(obsData[, mmVariables])
meanData <- attr(normData, "scaled:center")
sdData   <- attr(normData, "scaled:scale")

# ==============================================================================
# Step 4: Estimate PLS-SEM Model
#
# Run the PLS-SEM algorithm via the seminr package's simplePLS function.
# - maxIt = 300: maximum number of iterations for the PLS algorithm
# - stopCriterion = 7: convergence threshold (10^-7 for weight change)
# - inner_weights = path_weighting: use the path weighting scheme
#
# After estimation, attach the standardized data and raw data to the model
# object so that evaluation functions can access them.
# ==============================================================================

maxIt         <- 300
stopCriterion <- 7

plsModel <- seminr::simplePLS(
  obsData = normData,
  smMatrix = smMatrix,
  mmMatrix = mmMatrix,
  inner_weights = path_weighting,
  maxIt = maxIt,
  stopCriterion = stopCriterion,
  measurement_mode_scheme = measurement_mode_scheme
)

plsModel$data              <- normData
plsModel$rawdata           <- obsData
plsModel$measurement_model <- mm_model
plsModel$structural_model  <- sm_model

# Extract the measurement variable names for the endogenous construct BI
measurement_name <- plsModel$mmMatrix[plsModel$mmMatrix[, "construct"] == "BI", "measurement"]

# ==============================================================================
# Step 5: Evaluate Model with PLSpredict (10-fold Cross-Validation)
#
# Use K-fold cross-validation (K=10) to assess out-of-sample predictive power.
# For each fold, PLS and a naive linear model (LM) are both estimated, and
# three prediction error metrics are computed:
#   - RMSE (Root Mean Square Error)
#   - MAPE (Mean Absolute Percentage Error)
#   - MAD  (Mean Absolute Deviation)
#
# If PLS metrics are lower than LM metrics, the model has predictive power
# beyond a simple linear benchmark.
#
# Since BI's indicators are ILR-transformed (each pair corresponds to one
# original composition), we also aggregate metrics at the compositional block
# level: (BI1,BI2), (BI3,BI4), (BI5,BI6) each represent one composition.
# Block RMSE = sqrt(mean of squared per-column RMSE within the block).
# ==============================================================================

source("PLSpredict.R")

predictionMetrics <- validatePredict(normData, smMatrix, mmMatrix, noFolds = 10)

cat("PLS RMSE:\n"); print(predictionMetrics$PLSRMSE)
cat("LM  RMSE:\n"); print(predictionMetrics$LMRMSE)
cat("PLS MAPE:\n"); print(predictionMetrics$PLSMAPE)
cat("LM  MAPE:\n"); print(predictionMetrics$LMMAPE)
cat("PLS MAD:\n");  print(predictionMetrics$PLSMAD)
cat("LM  MAD:\n");  print(predictionMetrics$LMMAD)

PLSRMSE <- predictionMetrics$PLSRMSE
LMRMSE  <- predictionMetrics$LMRMSE
PLSMAPE <- predictionMetrics$PLSMAPE
LMMAPE  <- predictionMetrics$LMMAPE
PLSMAD  <- predictionMetrics$PLSMAD
LMMAD   <- predictionMetrics$LMMAD

# Define compositional blocks: each pair of ILR columns maps to one composition
blocks <- list(c("BI1", "BI2"), c("BI3", "BI4"), c("BI5", "BI6"))

calculate_block_metrics <- function(metrics_rmse, metrics_mape, metrics_mad, blocks) {
  results <- data.frame()
  for (block in blocks) {
    RMSE_block <- sqrt(rowSums(metrics_rmse[, block, drop = FALSE]^2) / length(block))
    MAPE_block <- rowMeans(metrics_mape[, block, drop = FALSE])
    MAD_block  <- rowMeans(metrics_mad[, block, drop = FALSE])
    results <- rbind(results, data.frame(
      block = paste0(block, collapse = ","),
      RMSE = mean(RMSE_block),
      MAPE = mean(MAPE_block),
      MAD  = mean(MAD_block)
    ))
  }
  return(results)
}

pls_block_metrics <- calculate_block_metrics(PLSRMSE, PLSMAPE, PLSMAD, blocks)
lm_block_metrics  <- calculate_block_metrics(LMRMSE, LMMAPE, LMMAD, blocks)

block_metrics_summary <- data.frame(
  block    = pls_block_metrics$block,
  RMSE_pls = pls_block_metrics$RMSE,
  RMSE_lm  = lm_block_metrics$RMSE,
  MAPE_pls = pls_block_metrics$MAPE,
  MAPE_lm  = lm_block_metrics$MAPE,
  MAD_pls  = pls_block_metrics$MAD,
  MAD_lm   = lm_block_metrics$MAD
)
print(block_metrics_summary)

# ==============================================================================
# Step 6: Evaluate Model Quality
#
# Measurement model evaluation:
#   - Cronbach's Alpha: internal consistency reliability for each construct.
#   - Composite Reliability (CR): similar to alpha but uses outer loadings.
#   - Compositional CR (COMP-CR): CR computed at the compositional block level
#     by grouping ILR indicator pairs back to their original compositions.
#     'comp_group' maps each indicator row to its parent composition.
#
# Structural model evaluation:
#   - R-squared: variance explained in endogenous constructs.
#   - Bootstrap (nboot=10000): significance testing of path coefficients via
#     resampling. Returns bootstrap means, SDs, and confidence intervals.
# ==============================================================================

source("evaluation.R")

# Map each ILR indicator to its parent composition for COMP-CR calculation.
# For example, Attitude1 & Attitude2 both map to composition "Attitude2",
# Attitude3 & Attitude4 map to "Attitude3", etc.
comp_group <- matrix(rep(
  c(multi_items("Attitude", 2:5),
    multi_items("BI", c(1, 2, 4)),
    multi_items("PBC", 3:5),
    "SNnew"),
  each = 2
))
rownames(comp_group) <- rownames(plsModel$outer_loadings)

cat("\n--- Cronbach's Alpha ---\n")
print(cronbachs_alpha(plsModel))

cat("\n--- Composite Reliability (CR) ---\n")
print(calculate_CR(plsModel))

cat("\n--- Compositional CR (COMP-CR) ---\n")
print(calculate_CR(plsModel, comp_group = comp_group))

cat("\n--- R-squared ---\n")
print(calculate_R2(plsModel))

cat("\n--- Bootstrap (nboot=10000) ---\n")
print(calculate_bootstrap(plsModel, nboot = 10000, seed = 123))
