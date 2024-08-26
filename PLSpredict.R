

#simplePLS 
#Description: This library contains the functions utilized to run the PLS-PM 
# algorithm.

#Function that estimates the PLS-PM Model
evaluate_simplePLS <- function(obsData,smMatrix, mmMatrix, maxIt=300, stopCriterion=7){
  
  #Create list of Measurements Variables
  mmVariables <- mmMatrix[,"measurement"]
  
  #Create list of Latent Variables
  ltVariables <- unique(c(smMatrix[,1],smMatrix[,2]))
  
  #Extract and Normalize the measurements for the model
  normData <- scale(obsData[,mmVariables],TRUE,TRUE)
  
  #Extract Mean and Standard Deviation of measurements for future prediction
  meanData <- attr(normData, "scaled:center")
  sdData <- attr(normData, "scaled:scale")
  
  #Create a matrix of outer_weights
  outer_weights <- matrix(data=0,
                          nrow=length(mmVariables),
                          ncol=length(ltVariables),
                          dimnames = list(mmVariables,ltVariables))
  
  #Initialize outer_weights matrix with value 1 for each relationship in the measurement model
  for (i in 1:length(ltVariables))  {
    outer_weights [mmMatrix[mmMatrix[,"construct"]==ltVariables[i],
                            "measurement"],
                   ltVariables[i]] =1
  }
  
  #Create a matrix of inner paths
  #? inner_paths => inner_weights?
  inner_paths <- matrix(data=0,
                        nrow=length(ltVariables),
                        ncol=length(ltVariables),
                        dimnames = list(ltVariables,ltVariables))
  
  #Iterative Process Starts here
  for (iterations in 0:maxIt)  {
    
    #Estimate Factor Scores from Outter Path
    #? fscores <- normData%*%outer_weights
    fscores <- normData[,mmVariables]%*%outer_weights
    
    #Standarize Factor Scores
    fscores <- scale(fscores,TRUE,TRUE)
    
    #Estimate inner paths (symmetric matrix)
    for (i in 1:nrow(smMatrix))  {
      inner_paths[smMatrix[i,"source"],
                  smMatrix[i,"target"]] = cov(fscores[,smMatrix[i,"source"]],
                                              fscores[,smMatrix[i,"target"]])
      #? next step necessary?
      inner_paths[smMatrix[i,"target"],
                  smMatrix[i,"source"]] = cov(fscores[,smMatrix[i,"source"]],
                                              fscores[,smMatrix[i,"target"]])
    }
    
    #Estimate Factor Scores from Inner Path
    fscores<-fscores%*%inner_paths
    
    #Standarize Factor Scores
    fscores <- scale(fscores,TRUE,TRUE)
    
    #Save last outer_weights
    last_outer_weights <- outer_weights
    
    #Update outer_weights
    for (i in 1:length(ltVariables))  {
      
      #If the measurement model is Formative
      if(mmMatrix[mmMatrix[,"construct"]==ltVariables[i],"type"][1]=="F"){
        outer_weights[mmMatrix[mmMatrix[,"construct"]==ltVariables[i], "measurement"], ltVariables[i]] = 
          solve(cor(normData[,mmMatrix[mmMatrix[,"constructt"]==ltVariables[i],"measurement"]])) %*%
          cor(normData[,mmMatrix[mmMatrix[,"construct"]==ltVariables[i],"measurement"]],
              fscores[,ltVariables[i]])
      }
      
      #If the measurement model is Reflective
      if(mmMatrix[mmMatrix[,"construct"]==ltVariables[i],"type"][1]=="R"){
        outer_weights[mmMatrix[mmMatrix[,"construct"]==ltVariables[i], "measurement"], ltVariables[i]] = 
          cov(normData[,mmMatrix[mmMatrix[,"construct"]==ltVariables[i],"measurement"]],fscores[,ltVariables[i]])
      }
    }
    
    #Estimate Factor Scores from Outer Weights
    fscores <- normData[,mmVariables]%*%outer_weights
    
    #Standarize outer_weights
    for (i in 1:length(ltVariables))  {
      outer_weights [mmMatrix[mmMatrix[,"construct"]==ltVariables[i], "measurement"], ltVariables[i]] =
        outer_weights [mmMatrix[mmMatrix[,"construct"]==ltVariables[i], "measurement"], ltVariables[i]] / sd(fscores[,ltVariables[i]])
    }
    
    #Verify the stop criteria
    weightDiff <- sum(abs(outer_weights-last_outer_weights))
    if (weightDiff <(10^(-(stopCriterion))))
      break
    
  } #Finish Iterative Process
  
  #Estimate Factor Scores from Outter Path
  fscores <- normData[,mmVariables]%*%outer_weights
  
  
  #Initialize Matrix of Path Coefficients
  path_coef <- matrix(data=0,
                      nrow=length(ltVariables),
                      ncol=length(ltVariables),
                      dimnames = list(ltVariables,ltVariables))
  
  #Identify which variables have incoming paths
  dependant<-unique(smMatrix[,"target"])
  
  #We calculate a linear regresion for each dependant variable
  for (i in 1:length(dependant))  {
    
    #Indentify the independant variables
    independant<-smMatrix[smMatrix[,"target"]==dependant[i],"source"]
    
    #Solve the sistem of equations
    results<- solve(cor(fscores[,independant, drop=FALSE]),
                    cor(fscores[,independant], fscores[,dependant[i]]))
    
    #Transform to a generic vector
    coefficients <- as.vector(results)
    if(!is.null(rownames(results)))
      names(coefficients)<-rownames(results)
    else
      names(coefficients)<-names(results)
    
    #Assign the Beta Values to the Path Coefficient Matrix
    for (j in 1:length(independant))  
      path_coef[independant[j],dependant[i]]=coefficients[independant[j]]
    
  }
  
  #Create a matrix of Outer Loadings
  outer_loadings <- matrix(data=0,
                           nrow=length(mmVariables),
                           ncol=length(ltVariables),
                           dimnames = list(mmVariables,ltVariables))
  
  
  #Calculate the Outer Loadings
  for (i in 1:length(ltVariables))  {
    outer_loadings [mmMatrix[mmMatrix[,"construct"]==ltVariables[i],
                             "measurement"],
                    ltVariables[i]] = cov(normData[,mmMatrix[mmMatrix[,"construct"]==ltVariables[i],"measurement"]],fscores[,ltVariables[i]])
    
  }
  
  #Calculate R Squared
  
  #Get smMatrix
  modelMatrix <- data.frame(smMatrix)
  
  #Get endogenous composites
  uniquetarget <- as.character(unique(modelMatrix$target)) 
  
  #Get composite scores
  valuesMatrix <- fscores
  
  #Calculate Linear Models
  lmmodels <- lapply(uniquetarget, function(x) {lm(as.formula(paste(x,"~ .", sep = "")), 
                                                   data = data.frame(valuesMatrix[,colnames(valuesMatrix) %in% 
                                                                                    c(x,as.character(modelMatrix$source[which(modelMatrix$target==x)]))]))})
  
  #Initialize matrix holder for Rsquared values
  rSquared <- matrix(,nrow=1,ncol=length(uniquetarget),byrow =TRUE,dimnames = list(1,uniquetarget))
  
  # Iterate and extract every R^2 value 
  for (i in 1:length(lmmodels)) {
    rSquared[,i] <- summary(lmmodels[[i]])$r.squared
  }
  
  
  #Prepare return Object
  plsModel <- list(meanData = meanData,
                   sdData = sdData,
                   smMatrix = smMatrix,
                   mmMatrix = mmMatrix,
                   ltVariables = ltVariables,
                   mmVariables = mmVariables,
                   outer_loadings = outer_loadings,
                   outer_weights = outer_weights,
                   path_coef = path_coef,
                   iterations = iterations,
                   weightDiff = weightDiff,
                   fscores = fscores,
                   rSquared = rSquared)
  
  class(plsModel) <- "plsModel"
  return(plsModel)
}


#Function that given a linear model, receives new data and uses the model to make predictions
predictlm <- function(model,newData){
  
  #Get Coefficients
  coef<-coefficients(model)
  
  #Initialize result vector with the intercept value
  result <-seq(coef["(Intercept)"],coef["(Intercept)"],length.out=nrow(newData))
  
  
  # Multiply the new values times their coefficient 
  for (i in 2:length(coef))  {
    result <- result+coef[i]*newData[names(coef)[i]]  
  }
  
  #Rename Column
  names(result)<-"result"
  
  #Return the prediction
  return (result)
}  



#PLSpredict
#Description: This library contains the functions utilized to run the PLS-PM 
# algorithm and its predictions.
# {
#   trainData = trainingData
#   testData = testingData
#   maxIt = 300
#   stopCriterion = 7
# }

#Function that receives a model and predicts measurements
PLSpredict <- function(trainData, testData, smMatrix, mmMatrix, maxIt=300, stopCriterion=7){
  
  #Call simplePLS function
  plsModel <- evaluate_simplePLS(trainData, smMatrix, mmMatrix, maxIt, stopCriterion)
  
  #Get results from model
  smMatrix <- plsModel$smMatrix
  mmMatrix <- plsModel$mmMatrix
  ltVariables <- plsModel$ltVariables
  mmVariables <- plsModel$mmVariables
  outer_weights <- plsModel$outer_weights
  outer_loadings <- plsModel$outer_loadings
  meanData<-plsModel$meanData
  sdData <- plsModel$sdData
  path_coef<-plsModel$path_coef
  
  #Create container for Exogenous Variables
  exVariables = NULL
  
  #Create container for Endogenous Variables
  enVariables = NULL
  
  #Identify Exogenous and Endogenous Variables
  exVariables <- unique(smMatrix[,1])
  pMeasurements <- NULL
  for (i in 1:length(exVariables)){
    pMeasurements <- c(pMeasurements,mmMatrix[mmMatrix[,"construct"]==exVariables[i],"measurement"])
  }
  enVariables <- unique(smMatrix[,2])
  resMeasurements <- NULL
  for (i in 1:length(enVariables)){
    resMeasurements <- c(resMeasurements, mmMatrix[mmMatrix[, "construct"] == enVariables[i],"measurement"])
  }
  enVariables <- setdiff(enVariables,exVariables)
  eMeasurements <- NULL
  for (i in 1:length(enVariables)){
    eMeasurements <- c(eMeasurements, mmMatrix[mmMatrix[, "construct"] == enVariables[i],"measurement"])
  }
  
  #Extract Measurements needed for Predictions
  normData <- testData[,pMeasurements]
  
  #Normalize data
  for (i in pMeasurements)
  {
    normData[,i] <-(normData[,i] - meanData[i])/sdData[i]
  }  
  
  #Convert dataset to matrix
  normData<-data.matrix(normData)
  
  #Add empty columns to normData for the estimated measurements
  for (i in 1:length(eMeasurements))
  {
    normData = cbind(normData, seq(0,0,length.out =nrow(normData)))
    colnames(normData)[length(colnames(normData))]=eMeasurements[i]
  }
  
  #Estimate Factor Scores from Outter Path
  fscores <- normData%*%outer_weights
  
  #Estimate Factor Scores from Inner Path and complete Matrix
  fscores <- fscores + fscores%*%path_coef
  
  #Predict Measurements with loadings
  predictedMeasurements<-fscores%*% t(outer_loadings)
  
  #Denormalize data
  for (i in mmVariables)
  {
    predictedMeasurements[,i]<-(predictedMeasurements[,i] * sdData[i])+meanData[i]
  }  
  
  #Calculating the residuals
  residuals <- testData[,eMeasurements] - predictedMeasurements[,eMeasurements]
  
  # # 输出预测值与实际值
  # cat("\nPredicted Measurements:\n")
  # print(predictedMeasurements)
  # cat("\nActual Measurements:\n")
  # print(testData[, eMeasurements])
  
  #Prepare return Object
  predictResults <- list(testData = testData[,eMeasurements],
                         predictedMeasurements = predictedMeasurements[,eMeasurements],
                         residuals = residuals,
                         compositeScores = fscores)
  
  class(predictResults) <- "predictResults"
  return(predictResults)
}



# validatePredict
# Description: This library contains the function utilized to generate perform a k-fold validation 
# and subsequent calculation of prediction metrics (RMSE, MAPE, MAD) for PLS & LM
# 
# {testData = normData
#   maxIt=300
#   stopCriterion=7
#   noFolds=10}

#Function for Generating and Evaluating Out-of-sample predictions using 10-fold cross-validation
validatePredict <- function(testData, smMatrix, mmMatrix, maxIt=300, stopCriterion=7, noFolds=10){
  
  #Randomly shuffle the data
  testData <- testData[sample(nrow(testData)),]
  
  #Create 10 equally size folds
  folds <- cut(seq(1, nrow(testData)), breaks=noFolds, labels=FALSE)
  
  #Identify variables to be tested
  uniqueTarget <- unique(smMatrix[,2])
  items <- NULL
  for (i in 1:length(uniqueTarget)){
    items <- c(items, mmMatrix[mmMatrix[, "construct"] == uniqueTarget[i], "measurement"])
  }
  uniqueSource <- unique(smMatrix[,1])
  sources <- NULL
  for (i in 1:length(uniqueSource)){
    sources <- c(sources, mmMatrix[mmMatrix[, "construct"] == uniqueSource[i], "measurement"])
  }
  lmtarget <- ifelse(length(intersect(uniqueTarget, uniqueSource)) == 0, uniqueTarget, setdiff(uniqueTarget, uniqueSource))  
  targets <- NULL
  for (i in 1:length(lmtarget)){
    targets <- c(targets, mmMatrix[mmMatrix[, "construct"] == lmtarget[i], "measurement"])
  }
  
  # Initialize matrices for prediction metrics
  PLSRMSE <- matrix(0, nrow=1, ncol=length(targets), byrow=TRUE, dimnames=list(1, targets))
  PLSSSE <- matrix(0, nrow=noFolds, ncol=length(targets), byrow=TRUE, dimnames=list(1:noFolds, targets))
  LMRMSE <- matrix(0, nrow=1, ncol=length(targets), byrow=TRUE, dimnames=list(1, targets))
  LMSSSE <- matrix(0, nrow=noFolds, ncol=length(targets), byrow=TRUE, dimnames=list(1:noFolds, targets))
  
  PLSSAPE <- matrix(0, nrow=noFolds, ncol=length(targets), byrow=TRUE, dimnames=list(1:noFolds, targets))
  PLSMAPE <- matrix(0, nrow=1, ncol=length(targets), byrow=TRUE, dimnames=list(1, targets))
  LMMAPE <- matrix(0, nrow=1, ncol=length(targets), byrow=TRUE, dimnames=list(1, targets))
  LMSAPE <- matrix(0, nrow=noFolds, ncol=length(targets), byrow=TRUE, dimnames=list(1:noFolds, targets))
  
  PLSSAD <- matrix(0, nrow=noFolds, ncol=length(targets), byrow=TRUE, dimnames=list(1:noFolds, targets))
  PLSMAD <- matrix(0, nrow=1, ncol=length(targets), byrow=TRUE, dimnames=list(1, targets))
  LMMAD <- matrix(0, nrow=1, ncol=length(targets), byrow=TRUE, dimnames=list(1, targets))
  LMSAD <- matrix(0, nrow=noFolds, ncol=length(targets), byrow=TRUE, dimnames=list(1:noFolds, targets))
  
  # Extract the target and non-target variables for Linear Model
  independentMatrix <- testData[,sources]
  dependentMatrix <- as.matrix(testData[,targets])
  
  #Perform 10 fold cross validation
  for(i in 1:noFolds){
    #Segment your data by fold using the which() function 
    testIndexes <- which(folds==i, arr.ind=TRUE)
    testingData <- testData[testIndexes, ]
    trainingData <- testData[-testIndexes, ]
    indepTestData <- independentMatrix[testIndexes, ]
    indepTrainData <- independentMatrix[-testIndexes, ]
    depTestData <- as.matrix(dependentMatrix[testIndexes, ])
    depTrainData <- as.matrix(dependentMatrix[-testIndexes, ])
    
    #PLS model
    testHolder <- PLSpredict(trainingData, testingData, smMatrix, mmMatrix, maxIt, stopCriterion)
    
    #Initialize PLS residuals and actuals holder matrices
    temptest <- as.matrix(testHolder$testData)
    tempresiduals <- as.matrix(testHolder$residuals)
    colnames(temptest) <- targets
    colnames(tempresiduals) <- targets
    PLSactuals <- as.matrix(temptest[, targets])
    PLSresiduals <- as.matrix(tempresiduals[, targets])
    
    # # 输出每次折叠的预测值与实际值
    # cat("\nFold", i, "PLS Predicted vs Actual:\n")
    # print(testHolder$predictedMeasurements)
    # print(testHolder$testData)

    #Initialize lm residuals and actuals holder matrices
    lmprediction <- matrix(0, nrow=nrow(depTestData), ncol=length(targets), byrow=TRUE, dimnames=list(1:nrow(depTestData), targets))
    lmresidual <- matrix(0, nrow=nrow(depTestData), ncol=length(targets), byrow=TRUE, dimnames=list(1:nrow(depTestData), targets))
    lmactual <- matrix(0, nrow=nrow(depTestData), ncol=length(targets), byrow=TRUE, dimnames=list(1:nrow(depTestData), targets))
    
    #LM Models
    for(l in 1:length(targets)){
      trainLM <- lm(depTrainData[, l] ~ ., as.data.frame(indepTrainData))
      lmprediction[, l] <- predict(trainLM, newdata=as.data.frame(indepTestData))
      lmresidual[, l] <- lmprediction[, l] - depTestData[, l]
      lmactual[, l] <- depTestData[, l]
    }

    # # 输出每次折叠的LM预测值与实际值
    # cat("\nFold", i, "LM Predicted vs Actual:\n")
    # print(lmprediction)
    # print(lmactual)

    #Iterate over no of targets
    for(j in 1:length(targets)){
      # Filter out near-zero actual values to avoid large MAPE values
      valid_indices_pls <- PLSactuals[, j] > 1e-6
      valid_PLSresiduals <- PLSresiduals[valid_indices_pls, j]
      valid_PLSactuals <- PLSactuals[valid_indices_pls, j]
      
      valid_indices_lm <- lmactual[, j] > 1e-6
      valid_lmresidual <- lmresidual[valid_indices_lm, j]
      valid_lmactual <- lmactual[valid_indices_lm, j]
      
      #Calculate SMSE
      PLSSSE[i, j] <- sum(valid_PLSresiduals^2)
      LMSSSE[i, j] <- sum(valid_lmresidual^2)
      #Calculate SAPE
      PLSSAPE[i, j] <- sum(abs(valid_PLSresiduals / valid_PLSactuals))
      LMSAPE[i, j] <- sum(abs(valid_lmresidual / valid_lmactual))
      #Calculate SAD
      PLSSAD[i, j] <- sum(abs(valid_PLSresiduals - mean(valid_PLSresiduals)))
      LMSAD[i, j] <- sum(abs(valid_lmresidual - mean(valid_lmresidual)))
    }
  }
  
  #Final calculations 
  denom <- noFolds * nrow(testData)  # Note: Changed from testingData to testData for correct denominator calculation
  for (k in 1:length(targets)) {
    LMRMSE[, k] <- sqrt((sum(LMSSSE[, k])) / denom)
    PLSRMSE[, k] <- sqrt((sum(PLSSSE[, k])) / denom)
    LMMAPE[, k] <- 100 * (sum(LMSAPE[, k]) / denom)
    PLSMAPE[, k] <- 100 * (sum(PLSSAPE[, k]) / denom)
    LMMAD[, k] <- sum(LMSAD[, k]) / denom
    PLSMAD[, k] <- sum(PLSSAD[, k]) / denom
  }
  
  validateResults <- list(
    PLSRMSE = PLSRMSE, 
    PLSMAPE = PLSMAPE,
    PLSMAD = PLSMAD,
    LMRMSE = LMRMSE,
    LMMAPE = LMMAPE,
    LMMAD = LMMAD
  )
  return(validateResults)
}
