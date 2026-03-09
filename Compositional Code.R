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
# Step 1: Data Loading and ILR Transformation
# ==============================================================================

data <- read_excel("TPB.xlsx", sheet = 1)
SN   <- data[, 31:32]
data <- data[, 1:30]

colnames(data) <- c(
  multi_items("Attitude", 1:12),
  multi_items("PBC", 1:9),
  multi_items("BI", 1:9)
)
data <- as.matrix(data)
data[data == 0] <- 0.001

becomp <- function(i) Closure(as.data.frame(data)[, (3 * i - 2):(3 * i)])
data <- do.call("cbind", apply(matrix(1:10), 1, becomp))

do.ilr <- function(i) ilr(data[, (3 * i - 2):(3 * i)])
data.ilr <- do.call("cbind", lapply(matrix(1:10), do.ilr))

colnames(data.ilr) <- c(
  multi_items("Attitude", 1:8),
  multi_items("PBC", 1:6),
  multi_items("BI", 1:6)
)
data.ilr <- cbind(data.ilr, SN)

head(data.ilr)
dim(data.ilr)

# ==============================================================================
# Step 2: PLS-SEM Model Estimation
# ==============================================================================

maxIt         <- 300
stopCriterion <- 7
obsData       <- data.ilr

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

normData <- standardize_safely(obsData[, mmVariables])
meanData <- attr(normData, "scaled:center")
sdData   <- attr(normData, "scaled:scale")

measurement_mode_scheme <- sapply(
  unique(c(sm_model[, 1], sm_model[, 2])),
  get_measure_mode, mmMatrix, USE.NAMES = TRUE
)

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

cat("Path Coefficients:\n")
print(plsModel$path_coef)

# ==============================================================================
# Step 3: Measurement Model Evaluation
# ==============================================================================

source("evaluation.R")

# Compositional loadings via ilr inverse transformation
write.csv(plsModel$outer_loadings, "real_outer_loadings.csv")
dd <- read.csv("real_outer_loadings-manual.csv")
dd1 <- data.frame()
for (i in 1:11) {
  dd1 <- rbind(dd1, ilr_i(matrix(as.numeric(dd[(2 * i - 1):(2 * i), 2]), ncol = 2)))
}
cat("\nCompositional Loadings:\n")
print(dd1)

# Compositional Composite Reliability
comp_group <- matrix(rep(
  c(multi_items("Attitude", 2:5),
    multi_items("BI", c(1, 2, 4)),
    multi_items("PBC", 3:5),
    "SNnew"),
  each = 2
))
rownames(comp_group) <- rownames(plsModel$outer_loadings)

cat("\n--- Compositional CR ---\n")
print(calculate_CR(plsModel, comp_group = comp_group))

cat("\n--- Cronbach's Alpha ---\n")
print(cronbachs_alpha(plsModel))

cat("\n--- Composite Reliability (CR) ---\n")
print(calculate_CR(plsModel))

# ==============================================================================
# Step 4: Structural Model Evaluation
# ==============================================================================

cat("\n--- Bootstrap ---\n")
print(calculate_bootstrap(plsModel, seed = 123))

cat("\n--- R-squared ---\n")
print(calculate_R2(plsModel))

cat("\n--- Q-squared ---\n")
print(calculate_Q2(plsModel, D = 10))

# ==============================================================================
# Step 5: Predictive Power Evaluation (10-fold CV)
# ==============================================================================

calculate_latent_scores <- function(outer_loadings, data) {
  latent_scores <- as.matrix(data) %*% ginv(t(as.matrix(outer_loadings)))
  colnames(latent_scores) <- colnames(outer_loadings)
  return(latent_scores)
}

calculate_predicted_endogenous_LV <- function(model, data, endogenous_LV) {
  latent_scores <- calculate_latent_scores(model$outer_loadings, data)
  latent_scores %*% model$path_coef[, endogenous_LV, drop = FALSE]
}

calculate_predicted_MV <- function(model, data, endogenous_LV) {
  measurement_name <- model$mmMatrix[model$mmMatrix[, "construct"] == endogenous_LV, "measurement"]
  latent_scores <- calculate_latent_scores(model$outer_loadings, data)

  predicted_MV <- matrix(0, nrow(latent_scores), length(measurement_name),
                         dimnames = list(NULL, measurement_name))
  for (i in 1:nrow(latent_scores)) {
    predicted_MV[i, ] <- latent_scores[i, endogenous_LV] *
      model$outer_loadings[measurement_name, endogenous_LV]
  }
  return(predicted_MV)
}

measurement_name <- plsModel$mmMatrix[plsModel$mmMatrix[, "construct"] == "BI", "measurement"]

# 10-fold cross-validation
folds <- createFolds(1:nrow(normData), k = 10)
n_folds <- length(folds)

RMSE_pls_LV <- numeric(n_folds)
MAE_pls_LV  <- numeric(n_folds)
RMSE_lm1_LV <- numeric(n_folds)
MAE_lm1_LV  <- numeric(n_folds)
RMSE_lm2_LV <- numeric(n_folds)
MAE_lm2_LV  <- numeric(n_folds)
RMSE_pls_MV <- matrix(0, n_folds, 6, dimnames = list(NULL, measurement_name))
MAE_pls_MV  <- matrix(0, n_folds, 6, dimnames = list(NULL, measurement_name))
RMSE_lm_MV  <- matrix(0, n_folds, 6, dimnames = list(NULL, measurement_name))
MAE_lm_MV   <- matrix(0, n_folds, 6, dimnames = list(NULL, measurement_name))

normData <- standardize_safely(obsData[, mmVariables])

for (i in seq_along(folds)) {
  test_indices  <- folds[[i]]
  train_indices <- setdiff(1:nrow(normData), test_indices)
  train_data <- plsModel$rawdata[train_indices, ]
  test_data  <- plsModel$rawdata[test_indices, ]

  plsModel_train <- seminr::simplePLS(
    obsData = train_data,
    smMatrix = smMatrix, mmMatrix = mmMatrix,
    inner_weights = path_weighting,
    maxIt = maxIt, stopCriterion = stopCriterion,
    measurement_mode_scheme = measurement_mode_scheme
  )

  # PLS: LV-level prediction
  estimated_test_eta <- calculate_predicted_endogenous_LV(plsModel_train, test_data, "BI")
  real_test_eta <- calculate_latent_scores(plsModel_train$outer_loadings, test_data)[, "BI"]

  RMSE_pls_LV[i] <- rmse(real_test_eta, estimated_test_eta)
  MAE_pls_LV[i]  <- mae(real_test_eta, estimated_test_eta)

  # LM1: LV-level regression (construct scores as predictors)
  train_latent_scores <- calculate_latent_scores(plsModel_train$outer_loadings, train_data)
  lm1 <- lm(BI ~ SNnew + Attitude + PBC, data = as.data.frame(train_latent_scores))

  test_latent_scores <- calculate_latent_scores(plsModel_train$outer_loadings, test_data)
  predicted_lm1 <- predict(lm1, newdata = as.data.frame(test_latent_scores))

  RMSE_lm1_LV[i] <- rmse(test_latent_scores[, "BI"], predicted_lm1)
  MAE_lm1_LV[i]  <- mae(test_latent_scores[, "BI"], predicted_lm1)

  # LM2: LV-level regression (MVs as predictors for LV)
  train_data_BI <- cbind(train_data, BI = train_latent_scores[, "BI"])
  test_data_BI  <- cbind(test_data,  BI = test_latent_scores[, "BI"])
  non_bi_vars   <- colnames(normData)[-(9:14)]
  lm2_formula   <- as.formula(paste("BI ~", paste(non_bi_vars, collapse = " + ")))
  lm2 <- lm(lm2_formula, data = train_data_BI)
  predicted_lm2 <- predict(lm2, newdata = test_data)

  real_test_eta <- calculate_latent_scores(plsModel_train$outer_loadings, test_data)[, "BI"]
  RMSE_lm2_LV[i] <- rmse(real_test_eta, predicted_lm2)
  MAE_lm2_LV[i]  <- mae(real_test_eta, predicted_lm2)

  # MV-level predictions
  estimated_MV <- calculate_predicted_MV(plsModel_train, test_data, "BI")
  actual_MV    <- test_data[, measurement_name]
  for (l in measurement_name) {
    RMSE_pls_MV[i, l] <- rmse(actual_MV[, l], estimated_MV[, l])
    MAE_pls_MV[i, l]  <- mae(actual_MV[, l], estimated_MV[, l])
  }

  # LM: MV-level regression
  train_data_X <- train_data[, !(colnames(train_data) %in% measurement_name)]
  test_data_X  <- test_data[, !(colnames(test_data) %in% measurement_name)]
  for (l in measurement_name) {
    train_data_l <- cbind(train_data_X, l = train_data[, l])
    test_data_l  <- cbind(test_data_X,  l = test_data[, l])
    trainLM <- lm(l ~ ., train_data_l)
    predicted_col <- predict(trainLM, newdata = test_data_l)
    RMSE_lm_MV[i, l] <- rmse(test_data[, l], predicted_col)
    MAE_lm_MV[i, l]  <- mae(test_data[, l], predicted_col)
  }
}

# --- MV-level Results ---
mat_MV <- rbind(
  colMeans(RMSE_pls_MV),
  colMeans(RMSE_lm_MV),
  colMeans(MAE_pls_MV),
  colMeans(MAE_lm_MV)
)

# Aggregate ilr pairs into compositional block metrics
mat_MV_final <- matrix(NA, nrow = 4, ncol = 3)
for (i in 1:2) {
  for (j in 1:3) {
    mat_MV_final[i, j] <- sqrt(0.5 * (mat_MV[i, 2 * j - 1]^2 + mat_MV[i, 2 * j]^2))
  }
}
for (i in 3:4) {
  for (j in 1:3) {
    mat_MV_final[i, j] <- 0.5 * (mat_MV[i, 2 * j - 1] + mat_MV[i, 2 * j])
  }
}
rownames(mat_MV_final) <- c("RMSE_PLS_MV", "RMSE_LM_MV", "MAE_PLS_MV", "MAE_LM_MV")
colnames(mat_MV_final) <- c("BI1", "BI2", "BI3")

cat("\n--- MV-level Predictions ---\n")
print(round(mat_MV_final, 3))

# --- LV-level Results ---
mat_LV <- cbind(
  apply(rbind(RMSE_pls_LV, RMSE_lm1_LV, RMSE_lm2_LV), 1, mean),
  apply(rbind(MAE_pls_LV, MAE_lm1_LV, MAE_lm2_LV), 1, mean)
)
colnames(mat_LV) <- c("RMSE_LV", "MAE_LV")
rownames(mat_LV) <- c("PLS", "LM1", "LM2")

cat("\n--- LV-level Predictions ---\n")
print(round(mat_LV, 3))
