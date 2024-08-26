rm(list = ls())
source('comp_function.R')
source('seminr.R')

library(readxl)
library(dplyr)
library(psych)
# library(seminr)
library(caret)
library(Metrics)
library(lmtest)
library(MASS)

######## Step 0: Data preprocessing ########
data = read_excel("data_sub.xlsx", sheet= 1)
SN = data[,31:32]
data = data[,1:30]

colnames(data) = c(multi_items("Attitude", 1:12),
                   multi_items("PBC", 1:9), 
                   multi_items("BI", 1:9)
                  )

data = as.matrix(data)
data[which(data == 0)] = 0.001

######## Step 1: Transform each compositional variable with ILR ########
# Define a function to apply compositional closure to sets of 3-part
becomp = function(i) {
  return(Closure(as.data.frame(data)[,(3*i-2):(3*i)]))
}

# Define a function to apply ilr transformation to sets of 3 columns
do.ilr = function(i) {
  return(ilr(data[,(3*i-2):(3*i)]))
}

# Apply the ilr transformation to the data
data = do.call("cbind", apply(matrix(1:10), 1, becomp))
data.ilr = do.call("cbind", lapply(matrix(1:10), do.ilr))

# Rename columns for the ilr-transformed data
colnames(data.ilr) = c(multi_items("Attitude", 1:8), 
                       multi_items("PBC", 1:6), 
                       multi_items("BI", 1:6)
                      )

# Combine the ilr-transformed data with the SN columns
data.ilr = cbind(data.ilr, SN)

head(data.ilr)
dim(data.ilr)
obsData <- data.ilr


######## Step 2: Specify the ILR-transformed indicators to relate to their constructs ########
######## as well as the measurement model and structural model ######## 
# Define the measurement model (outer model) with constructs and their corresponding indicators
mm_model <- constructs(
  composite("Attitude", multi_items("Attitude", 1:8), weights = mode_A),
  composite("BI", multi_items("BI", 1:6), weights = mode_A),
  composite("PBC", multi_items("PBC", 1:6), weights = mode_A),
  composite("SNnew", multi_items("SNnew", 1:2), weights = mode_A)
)

# Define the structural model (inner model) with paths between constructs
sm_model <- relationships(
  paths(from = "SNnew", to = "Attitude"),
  paths(from = "SNnew", to = "PBC"),
  paths(from = "SNnew", to = "BI"), 
  paths(from = "Attitude", to = "BI"),
  paths(from = "PBC", to = "BI")
)

# Create the structural model matrix
smMatrix <- sm_model
measurement_model <- mm_model

# Extract constructs and measurements for the measurement model
recognized_constructs <- c("composite", "reflective", "higher_order_composite")
construct_measurements <- measurement_model[names(measurement_model) %in% recognized_constructs]
mmMatrix <- matrix(unlist(construct_measurements),
                   ncol = 3, byrow = TRUE,
                   dimnames = list(NULL, c("construct", "measurement", "type"))
)

# Get the construct names from the structural model matrix
constructs <- construct_names(smMatrix)
# Extract measurement variables for the constructs
mmVariables <- mmMatrix[mmMatrix[,"construct"] %in% constructs, "measurement"]
# Determine the measurement mode scheme for each construct
measurement_mode_scheme <- sapply(unique(c(sm_model[,1], sm_model[,2])), get_measure_mode, mmMatrix, USE.NAMES = TRUE)


######## Step 3: Standardize the data ########
# Standardize the data
normData <- standardize_safely(obsData[, mmVariables])
meanData <- attr(normData, "scaled:center")
sdData <- attr(normData, "scaled:scale")

######## Step 4: Estimate using the standard PLS-SEM algorithm ########
# Set the maximum number of iterations and stop criterion for the PLS algorithm
maxIt <- 300
stopCriterion <- 7
# Fit the PLS-SEM model using seminr packages
plsModel = seminr::simplePLS(obsData = normData,
                             smMatrix = smMatrix,
                             mmMatrix = mmMatrix,
                             inner_weights = path_weighting,
                             maxIt = maxIt,
                             stopCriterion = stopCriterion,
                             measurement_mode_scheme = measurement_mode_scheme)

# Add standardized and raw data to the plsModel object
plsModel$data = normData
plsModel$rawdata = obsData
plsModel$measurement_model = mm_model
plsModel$structural_model = sm_model

# Extract measurement variables name
measurement_name <- plsModel$mmMatrix[plsModel$mmMatrix[,"construct"] == "BI","measurement"]

######## Step 5: Evaluate the model with PLSpredict ######## 
source('PLSpredict.R')
# Use Prediction Metrics in CVPAT Method
predictionMetrics <- validatePredict(normData, smMatrix, mmMatrix, noFolds=10)
predictionMetrics$PLSRMSE
predictionMetrics$LMRMSE
predictionMetrics$PLSMAPE
predictionMetrics$LMMAPE
predictionMetrics$PLSMAD
predictionMetrics$LMMAD

# Load prediction metrics from the CVPAT results
PLSRMSE <- predictionMetrics$PLSRMSE
LMRMSE <- predictionMetrics$LMRMSE
PLSMAPE <- predictionMetrics$PLSMAPE
LMMAPE <- predictionMetrics$LMMAPE
PLSMAD <- predictionMetrics$PLSMAD
LMMAD <- predictionMetrics$LMMAD

# Define blocks for the compositional data
blocks <- list(
  c("BI1", "BI2"),
  c("BI3", "BI4"),
  c("BI5", "BI6")
)

# Function to calculate block metrics (RMSE, MAPE, MAD)
calculate_block_metrics <- function(metrics_rmse, metrics_mape, metrics_mad, blocks) {
  results <- data.frame()
  
  # Iterate through each block to calculate metrics
  for (block in blocks) {
    rmse_cols <- block
    mape_cols <- block
    mad_cols <- block
    
    # Calculate RMSE, MAPE, and MAD for the block
    RMSE_block <- sqrt(rowSums(metrics_rmse[, rmse_cols, drop = FALSE]^2) / length(rmse_cols))
    MAPE_block <- rowMeans(metrics_mape[, mape_cols, drop = FALSE])
    MAD_block <- rowMeans(metrics_mad[, mad_cols, drop = FALSE])
    
    # Store the calculated metrics in the results dataframe
    results <- rbind(results, data.frame(
      block = paste0(block, collapse = ","),
      RMSE = mean(RMSE_block),
      MAPE = mean(MAPE_block),
      MAD = mean(MAD_block)
    ))
  }
  
  return(results)
}

# Calculate block metrics for PLS model
pls_block_metrics <- calculate_block_metrics(PLSRMSE, PLSMAPE, PLSMAD, blocks)

# Calculate block metrics for LM model
lm_block_metrics <- calculate_block_metrics(LMRMSE, LMMAPE, LMMAD, blocks)

# Combine the results for PLS and LM models into a summary dataframe
block_metrics_summary <- data.frame(
  block = pls_block_metrics$block,
  RMSE_pls = pls_block_metrics$RMSE,
  RMSE_lm = lm_block_metrics$RMSE,
  MAPE_pls = pls_block_metrics$MAPE,
  MAPE_lm = lm_block_metrics$MAPE,
  MAD_pls = pls_block_metrics$MAD,
  MAD_lm = lm_block_metrics$MAD
)

# Print the summary of block metrics
print(block_metrics_summary)

######## Step 6: Evaluate the model using measures such as COMP-CR, Cronbach's alpha, and etc. ######## 
# Source the evaluation functions from evaluation.R
source('evaluation.R')

# Define the compositional group for the evaluation metrics
# 'comp_group' matrix contains the groupings of measurement items for compositional blocks
comp_group = matrix(rep(c(multi_items('Attitude', 2:5), 
                          multi_items('BI', c(1, 2, 4)), 
                          multi_items('PBC', 3:5), 
                          'SNnew'), each = 2))
rownames(comp_group) = rownames(plsModel$outer_loadings)

# Calculate and print Cronbach's alpha
print('cronbachs_alpha')
cronbachs_alpha(plsModel)

# Calculate and print COMP-CR
print('CR')
calculate_CR(plsModel)
calculate_CR(plsModel, comp_group = comp_group) # With compositional group

# Calculate and print R-squared (R2)
print('R2')
calculate_R2(plsModel)

# Calculate and print Bootstrap confidence intervals with 10000 resamples
print('bootstrap')
calculate_bootstrap(plsModel, nboot = 10000, seed = 123)

