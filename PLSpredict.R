# ==============================================================================
# PLSpredict.R - PLS-PM estimation, prediction, and k-fold cross-validation
# ==============================================================================

# --- PLS-PM Model Estimation ---

evaluate_simplePLS <- function(obsData, smMatrix, mmMatrix, maxIt = 300, stopCriterion = 7) {
  mmVariables <- mmMatrix[, "measurement"]
  ltVariables <- unique(c(smMatrix[, 1], smMatrix[, 2]))

  normData <- scale(obsData[, mmVariables], TRUE, TRUE)
  meanData <- attr(normData, "scaled:center")
  sdData   <- attr(normData, "scaled:scale")

  outer_weights <- matrix(0, nrow = length(mmVariables), ncol = length(ltVariables),
                          dimnames = list(mmVariables, ltVariables))

  for (i in seq_along(ltVariables)) {
    outer_weights[mmMatrix[mmMatrix[, "construct"] == ltVariables[i], "measurement"],
                  ltVariables[i]] <- 1
  }

  inner_paths <- matrix(0, nrow = length(ltVariables), ncol = length(ltVariables),
                        dimnames = list(ltVariables, ltVariables))

  # Iterative PLS algorithm
  for (iterations in 0:maxIt) {
    fscores <- scale(normData[, mmVariables] %*% outer_weights, TRUE, TRUE)

    for (i in 1:nrow(smMatrix)) {
      src <- smMatrix[i, "source"]
      tgt <- smMatrix[i, "target"]
      inner_paths[src, tgt] <- cov(fscores[, src], fscores[, tgt])
      inner_paths[tgt, src] <- inner_paths[src, tgt]
    }

    fscores <- scale(fscores %*% inner_paths, TRUE, TRUE)

    last_outer_weights <- outer_weights

    for (i in seq_along(ltVariables)) {
      items <- mmMatrix[mmMatrix[, "construct"] == ltVariables[i], "measurement"]
      type  <- mmMatrix[mmMatrix[, "construct"] == ltVariables[i], "type"][1]

      if (type == "F") {
        outer_weights[items, ltVariables[i]] <-
          solve(cor(normData[, items])) %*% cor(normData[, items], fscores[, ltVariables[i]])
      }
      if (type == "R") {
        outer_weights[items, ltVariables[i]] <-
          cov(normData[, items], fscores[, ltVariables[i]])
      }
    }

    fscores <- normData[, mmVariables] %*% outer_weights

    for (i in seq_along(ltVariables)) {
      items <- mmMatrix[mmMatrix[, "construct"] == ltVariables[i], "measurement"]
      outer_weights[items, ltVariables[i]] <-
        outer_weights[items, ltVariables[i]] / sd(fscores[, ltVariables[i]])
    }

    weightDiff <- sum(abs(outer_weights - last_outer_weights))
    if (weightDiff < 10^(-stopCriterion)) break
  }

  fscores <- normData[, mmVariables] %*% outer_weights

  # Path coefficients
  path_coef <- matrix(0, nrow = length(ltVariables), ncol = length(ltVariables),
                      dimnames = list(ltVariables, ltVariables))
  dependant <- unique(smMatrix[, "target"])

  for (i in seq_along(dependant)) {
    independant <- smMatrix[smMatrix[, "target"] == dependant[i], "source"]
    results <- solve(cor(fscores[, independant, drop = FALSE]),
                     cor(fscores[, independant], fscores[, dependant[i]]))
    coefficients <- as.vector(results)
    names(coefficients) <- if (!is.null(rownames(results))) rownames(results) else names(results)

    for (j in seq_along(independant)) {
      path_coef[independant[j], dependant[i]] <- coefficients[independant[j]]
    }
  }

  # Outer loadings
  outer_loadings <- matrix(0, nrow = length(mmVariables), ncol = length(ltVariables),
                           dimnames = list(mmVariables, ltVariables))
  for (i in seq_along(ltVariables)) {
    items <- mmMatrix[mmMatrix[, "construct"] == ltVariables[i], "measurement"]
    outer_loadings[items, ltVariables[i]] <-
      cov(normData[, items], fscores[, ltVariables[i]])
  }

  # R-squared via linear models
  modelMatrix <- data.frame(smMatrix)
  uniquetarget <- as.character(unique(modelMatrix$target))
  valuesMatrix <- fscores

  lmmodels <- lapply(uniquetarget, function(x) {
    predictors <- as.character(modelMatrix$source[modelMatrix$target == x])
    lm(as.formula(paste(x, "~ .", sep = "")),
       data = data.frame(valuesMatrix[, colnames(valuesMatrix) %in% c(x, predictors)]))
  })

  rSquared <- matrix(, nrow = 1, ncol = length(uniquetarget),
                     byrow = TRUE, dimnames = list(1, uniquetarget))
  for (i in seq_along(lmmodels)) {
    rSquared[, i] <- summary(lmmodels[[i]])$r.squared
  }

  plsModel <- list(
    meanData = meanData, sdData = sdData,
    smMatrix = smMatrix, mmMatrix = mmMatrix,
    ltVariables = ltVariables, mmVariables = mmVariables,
    outer_loadings = outer_loadings, outer_weights = outer_weights,
    path_coef = path_coef, iterations = iterations,
    weightDiff = weightDiff, fscores = fscores, rSquared = rSquared
  )
  class(plsModel) <- "plsModel"
  return(plsModel)
}

# --- LM Prediction Helper ---

predictlm <- function(model, newData) {
  coef <- coefficients(model)
  result <- rep(coef["(Intercept)"], nrow(newData))
  for (i in 2:length(coef)) {
    result <- result + coef[i] * newData[[names(coef)[i]]]
  }
  names(result) <- "result"
  return(result)
}

# --- PLS Prediction ---

PLSpredict <- function(trainData, testData, smMatrix, mmMatrix, maxIt = 300, stopCriterion = 7) {
  plsModel <- evaluate_simplePLS(trainData, smMatrix, mmMatrix, maxIt, stopCriterion)

  ltVariables    <- plsModel$ltVariables
  mmVariables    <- plsModel$mmVariables
  outer_weights  <- plsModel$outer_weights
  outer_loadings <- plsModel$outer_loadings
  meanData       <- plsModel$meanData
  sdData         <- plsModel$sdData
  path_coef      <- plsModel$path_coef

  # Identify exogenous/endogenous variables and their measurements
  exVariables <- unique(smMatrix[, 1])
  pMeasurements <- unlist(lapply(exVariables, function(v) {
    mmMatrix[mmMatrix[, "construct"] == v, "measurement"]
  }))

  enVariables <- setdiff(unique(smMatrix[, 2]), exVariables)
  eMeasurements <- unlist(lapply(enVariables, function(v) {
    mmMatrix[mmMatrix[, "construct"] == v, "measurement"]
  }))

  # Normalize test data using training parameters
  normData <- testData[, pMeasurements]
  for (i in pMeasurements) {
    normData[, i] <- (normData[, i] - meanData[i]) / sdData[i]
  }
  normData <- data.matrix(normData)

  # Add placeholder columns for endogenous measurements
  for (i in seq_along(eMeasurements)) {
    normData <- cbind(normData, rep(0, nrow(normData)))
    colnames(normData)[ncol(normData)] <- eMeasurements[i]
  }

  fscores <- normData %*% outer_weights
  fscores <- fscores + fscores %*% path_coef

  predictedMeasurements <- fscores %*% t(outer_loadings)

  # Denormalize predictions
  for (i in mmVariables) {
    predictedMeasurements[, i] <- predictedMeasurements[, i] * sdData[i] + meanData[i]
  }

  residuals <- testData[, eMeasurements] - predictedMeasurements[, eMeasurements]

  predictResults <- list(
    testData = testData[, eMeasurements],
    predictedMeasurements = predictedMeasurements[, eMeasurements],
    residuals = residuals,
    compositeScores = fscores
  )
  class(predictResults) <- "predictResults"
  return(predictResults)
}

# --- K-Fold Cross-Validation ---

validatePredict <- function(testData, smMatrix, mmMatrix, maxIt = 300, stopCriterion = 7, noFolds = 10) {
  testData <- testData[sample(nrow(testData)), ]
  folds <- cut(seq_len(nrow(testData)), breaks = noFolds, labels = FALSE)

  # Identify target measurements
  uniqueTarget <- unique(smMatrix[, 2])
  items <- unlist(lapply(uniqueTarget, function(v) {
    mmMatrix[mmMatrix[, "construct"] == v, "measurement"]
  }))

  uniqueSource <- unique(smMatrix[, 1])
  sources <- unlist(lapply(uniqueSource, function(v) {
    mmMatrix[mmMatrix[, "construct"] == v, "measurement"]
  }))

  lmtarget <- setdiff(uniqueTarget, uniqueSource)
  if (length(lmtarget) == 0) lmtarget <- uniqueTarget
  targets <- unlist(lapply(lmtarget, function(v) {
    mmMatrix[mmMatrix[, "construct"] == v, "measurement"]
  }))

  # Initialize metric accumulators
  init_mat <- function(nr, nc, nms) matrix(0, nrow = nr, ncol = nc, dimnames = list(1:nr, nms))

  PLSSSE  <- init_mat(noFolds, length(targets), targets)
  LMSSSE  <- init_mat(noFolds, length(targets), targets)
  PLSSAPE <- init_mat(noFolds, length(targets), targets)
  LMSAPE  <- init_mat(noFolds, length(targets), targets)
  PLSSAD  <- init_mat(noFolds, length(targets), targets)
  LMSAD   <- init_mat(noFolds, length(targets), targets)

  independentMatrix <- testData[, sources]
  dependentMatrix   <- as.matrix(testData[, targets])

  for (i in 1:noFolds) {
    testIndexes  <- which(folds == i)
    testingData  <- testData[testIndexes, ]
    trainingData <- testData[-testIndexes, ]

    indepTestData  <- independentMatrix[testIndexes, ]
    indepTrainData <- independentMatrix[-testIndexes, ]
    depTestData    <- as.matrix(dependentMatrix[testIndexes, ])
    depTrainData   <- as.matrix(dependentMatrix[-testIndexes, ])

    # PLS prediction
    testHolder   <- PLSpredict(trainingData, testingData, smMatrix, mmMatrix, maxIt, stopCriterion)
    PLSactuals   <- as.matrix(testHolder$testData)
    PLSresiduals <- as.matrix(testHolder$residuals)
    colnames(PLSactuals)   <- targets
    colnames(PLSresiduals) <- targets

    # LM prediction
    lmprediction <- matrix(0, nrow = nrow(depTestData), ncol = length(targets),
                           dimnames = list(NULL, targets))
    lmresidual   <- matrix(0, nrow = nrow(depTestData), ncol = length(targets),
                           dimnames = list(NULL, targets))
    lmactual     <- matrix(0, nrow = nrow(depTestData), ncol = length(targets),
                           dimnames = list(NULL, targets))

    for (l in seq_along(targets)) {
      trainLM <- lm(depTrainData[, l] ~ ., as.data.frame(indepTrainData))
      lmprediction[, l] <- predict(trainLM, newdata = as.data.frame(indepTestData))
      lmresidual[, l]   <- lmprediction[, l] - depTestData[, l]
      lmactual[, l]     <- depTestData[, l]
    }

    # Accumulate per-fold metrics
    for (j in seq_along(targets)) {
      valid_pls <- PLSactuals[, j] > 1e-6
      valid_lm  <- lmactual[, j] > 1e-6

      PLSSSE[i, j]  <- sum(PLSresiduals[valid_pls, j]^2)
      LMSSSE[i, j]  <- sum(lmresidual[valid_lm, j]^2)
      PLSSAPE[i, j] <- sum(abs(PLSresiduals[valid_pls, j] / PLSactuals[valid_pls, j]))
      LMSAPE[i, j]  <- sum(abs(lmresidual[valid_lm, j] / lmactual[valid_lm, j]))
      PLSSAD[i, j]  <- sum(abs(PLSresiduals[valid_pls, j] - mean(PLSresiduals[valid_pls, j])))
      LMSAD[i, j]   <- sum(abs(lmresidual[valid_lm, j] - mean(lmresidual[valid_lm, j])))
    }
  }

  # Final metric calculations
  denom <- noFolds * nrow(testData)
  PLSRMSE <- matrix(sqrt(colSums(PLSSSE) / denom), nrow = 1, dimnames = list(1, targets))
  LMRMSE  <- matrix(sqrt(colSums(LMSSSE) / denom),  nrow = 1, dimnames = list(1, targets))
  PLSMAPE <- matrix(100 * colSums(PLSSAPE) / denom,  nrow = 1, dimnames = list(1, targets))
  LMMAPE  <- matrix(100 * colSums(LMSAPE) / denom,   nrow = 1, dimnames = list(1, targets))
  PLSMAD  <- matrix(colSums(PLSSAD) / denom,         nrow = 1, dimnames = list(1, targets))
  LMMAD   <- matrix(colSums(LMSAD) / denom,          nrow = 1, dimnames = list(1, targets))

  list(
    PLSRMSE = PLSRMSE, PLSMAPE = PLSMAPE, PLSMAD = PLSMAD,
    LMRMSE = LMRMSE,   LMMAPE = LMMAPE,   LMMAD = LMMAD
  )
}
