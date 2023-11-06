constructs <- function(...) {
  return_list <- list(...)
  names(return_list) <- lapply(return_list, function(x) class(x)[[3]])
  class(return_list) <- c(class(return_list), "measurement_model", "seminr_model")
  return(return_list)
}


reflective <- function(construct_name, item_names) {
  construct_names <- rep(construct_name, length(item_names))
  construct <- c(rbind(construct_names, item_names, "C"))
  class(construct) <- c(class(construct), c("construct", "reflective"))
  return(construct)
}

composite <- function(construct_name, item_names, weights = correlation_weights) {
  if (identical(weights, correlation_weights)) {
    composite_type = "A"
  } else if (identical(weights, regression_weights)) {
    composite_type = "B"
  } else {
    stop("Composites must be defined as mode A (correlation weights) or B (regression weights)")
  }
  construct <- c(rbind(construct_name, item_names, composite_type))
  class(construct) <- c(class(construct), c("construct", "composite"))
  return(construct)
}

mode_A  <- function(mmMatrix, i, normData, construct_scores) {
  return(stats::cov(normData[,mmMatrix[mmMatrix[,"construct"]==i,"measurement"]],construct_scores[,i]))
}
#' @export
correlation_weights <- mode_A

mode_B <- function(mmMatrix, i,normData, construct_scores) {
  return(solve(stats::cor(normData[,mmMatrix[mmMatrix[,"construct"]==i,"measurement"]])) %*%
           stats::cor(normData[,mmMatrix[mmMatrix[,"construct"]==i,"measurement"]],
                      construct_scores[,i]))
}
#' @export
regression_weights <- mode_B

# measurement model
constructs <- function(...) {
  return_list <- list(...)
  names(return_list) <- lapply(return_list, function(x) class(x)[[3]])
  # class(return_list) <- c(class(return_list), "measurement_model", "seminr_model")
  return(return_list)
}

reflective <- function(construct_name, item_names) {
  construct_names <- rep(construct_name, length(item_names))
  construct <- c(rbind(construct_names, item_names, "C"))
  class(construct) <- c(class(construct), c("construct", "reflective"))
  return(construct)
}

composite <- function(construct_name, item_names, weights = correlation_weights) {
  if (identical(weights, correlation_weights)) {
    composite_type = "A"
  } else if (identical(weights, regression_weights)) {
    composite_type = "B"
  } else {
    stop("Composites must be defined as mode A (correlation weights) or B (regression weights)")
  }
  construct <- c(rbind(construct_name, item_names, composite_type))
  class(construct) <- c(class(construct), c("construct", "composite"))
  return(construct)
}

multi_items <- function(item_name, item_numbers, ...) {
  affix <- as.data.frame(list(...))
  paste(affix$prefix, item_name, affix$mid, item_numbers, affix$suffix, sep = "")
}

single_item <- function(item) {
  return(item)
}

mode_A  <- function(mmMatrix, i, normData, construct_scores) {
  return(stats::cov(normData[,mmMatrix[mmMatrix[,"construct"]==i,"measurement"]],construct_scores[,i]))
}

correlation_weights <- mode_A

mode_B <- function(mmMatrix, i,normData, construct_scores) {
  return(solve(stats::cor(normData[,mmMatrix[mmMatrix[,"construct"]==i,"measurement"]])) %*%
           stats::cor(normData[,mmMatrix[mmMatrix[,"construct"]==i,"measurement"]],
                      construct_scores[,i]))
}

regression_weights <- mode_B

# structural model 
relationships <- function(...) {
  smMatrix <- matrix(c(...), ncol = 2, byrow = TRUE,
                     dimnames = list(NULL, c("source", "target")))
  # class(smMatrix) <- c(class(smMatrix), "structural_model", "seminr_model")
  return(smMatrix)
}

paths <- function(from, to) {
  return(as.vector(t(as.matrix(expand.grid(from, to)))))
}

mmMatrix_per_construct <- function(construct, mmMatrix) {
  constructmatrix <- mmMatrix[mmMatrix[,"construct"]==construct,c("construct","measurement","type")]
  # If single item construct
  if (class(constructmatrix)[1] != "matrix") {
    constructmatrix = t(as.matrix(constructmatrix))
  }
  return(constructmatrix)
}

# rho[i, 1] <- compute_construct_rhoA(outer_weights, mmMatrix, construct = i, obsData)
compute_construct_rhoA <- function(weights, mmMatrix, construct, obsData) {
  # get the weights for the construct
  w <- as.matrix(weights[mmMatrix[mmMatrix[, "construct"]==construct, "measurement"], construct])
  
  # Get empirical covariance matrix of lv indicators (S)
  indicators <- scale(obsData[,mmMatrix[mmMatrix[,"construct"]==construct, "measurement"]], TRUE, TRUE)
  S <- stats::cov(indicators, indicators)
  diag(S) <- 0
  
  # Get AA matrix without diagonal
  AAnondiag <- w %*% t(w)
  diag(AAnondiag) <- 0
  
  # Calculate rhoA
  return((t(w) %*% w)^2 * ((t(w) %*% (S) %*% w)/(t(w) %*% AAnondiag %*% w)))
}

all_reflective <- function(mmMatrix, constructs) {
  unique(mmMatrix[mmMatrix[, "type"]=="C", "construct"])
}

interaction_term <- function(iv, moderator, method=product_indicator, weights = mode_A) {
  intxn <- method(iv, moderator, weights)
  class(intxn) <- class(method())
  
  return(intxn)
}

construct_names <- function(smMatrix) {
  unique(c(smMatrix[,1], smMatrix[,2]))
}

all_endogenous <- function(smMatrix) {
  unique(smMatrix[, "target"])
}


constructs(
  reflective("Image",        multi_items("IMAG", 1:5)),
  reflective("Expectation",  multi_items("CUEX", 1:3)),
  reflective("Quality",      multi_items("PERQ", 1:7)),
  reflective("Value",        multi_items("PERV", 1:2)),
  reflective("Satisfaction", multi_items("CUSA", 1:3)),
  reflective("Complaints",   single_item("CUSCO")),
  reflective("Loyalty",      multi_items("CUSL", 1:3))
)
constructs(
  composite("Image",        multi_items("IMAG", 1:5), weights = correlation_weights),
  composite("Expectation",  multi_items("CUEX", 1:3), weights = mode_A),
  composite("Quality",      multi_items("PERQ", 1:7), weights = regression_weights),
  composite("Value",        multi_items("PERV", 1:2), weights = mode_B)
)

# Distinguish and mix composite measurement (used in PLS-PM)
# or reflective (common-factor) measurement (used in CBSEM, CFA, and PLSc)
# - We will first use composites in PLS-PM analysis
# - Later we will convert the omposites into reflectives for CFA/CBSEM (step 3)
# measurements <- constructs(
#   composite("Image",        multi_items("IMAG", 1:5)),
#   composite("Expectation",  multi_items("CUEX", 1:3)),
#   composite("Value",        multi_items("PERV", 1:2)),
#   composite("Satisfaction", multi_items("CUSA", 1:3)),
#   interaction_term(iv = "Image", moderator = "Expectation")
# )

# Quickly create multiple paths "from" and "to" sets of constructs  
structure <- relationships(
  paths(from = c("Image", "Expectation", "Image*Expectation"), to = "Value"),
  paths(from = "Value", to = "Satisfaction")
)

standardize_safely <- function(x) {
  # NOTE: we could return zeros for columns with zero variance:
  # apply(x, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
  res <- scale(x, TRUE, TRUE)
  if (any(is.nan(res))) stop("zero variance items cannot be scaled")
  res
}

#' Inner weighting scheme functions to estimate inner paths matrix
path_factorial <- function(smMatrix,construct_scores, dependant, paths_matrix) {
  inner_paths <- stats::cor(construct_scores,construct_scores) * (paths_matrix + t(paths_matrix))
  return(inner_paths)
}


path_weighting <- function(smMatrix, construct_scores, dependant, paths_matrix) {
  # correlations for outgoing paths
  inner_paths <- stats::cor(construct_scores,construct_scores) * t(paths_matrix)
  
  #Regression betas for the incoming paths
  #Iterate and regress the incoming paths
  for (i in 1:length(dependant))  {
    #Indentify the independant variables
    independant<-smMatrix[smMatrix[,"target"]==dependant[i],"source"]
    
    #Solve the system of equations
    inner_paths[independant,dependant[i]] = solve(t(construct_scores[,independant]) %*% construct_scores[,independant], t(construct_scores[,independant]) %*% construct_scores[,dependant[i]])
  }
  return(inner_paths)
}

# function to get measurement mode of a construct (first item) as a function
get_measure_mode <- function(construct,mmMatrix) {
  if((mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "A")
     |(mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "C")
     |(mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "HOCA")) {
    return(mode_A)
  } else if((mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "B")
            |(mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "HOCB")) {
    return(mode_B)
  } else if((mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "UNIT") ) {
    return(unit_weights)
  }
  # ifelse((mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "A")
  #        |(mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "C")
  #        |(mmMatrix[mmMatrix[,"construct"]==construct,"type"][1] == "HOCA"), return(mode_A), return(mode_B))
}

standardize_outer_weights <- function(normData, mmVariables, outer_weights) {
  # Standardize the outer weights
  std_devs <- attr(scale((normData[,mmVariables]%*%outer_weights), center = FALSE),"scaled:scale")
  # divide matrix by std_devs and return
  return(t(t(outer_weights) / std_devs))
}

calculate_loadings <- function(weights_matrix,construct_scores, normData) {
  return(as.matrix(stats::cov(normData,construct_scores) * weights_matrix))
}

# Function to adjust for the interaction
# Adjustment of the SD of the interaction term as per Henseler, J., & Chin, W. W. (2010),
# A comparison of approaches for the analysis of interaction effects between latent variables
# using partial least squares path modeling. Structural Equation Modeling, 17(1), 82-109. https://doi.org/10.1080/10705510903439003
adjust_interaction <- function(constructs, mmMatrix, outer_loadings, construct_scores, obsData){
  for(construct in constructs) {
    adjustment <- 0
    denom <- 0
    if(grepl("\\*", construct)) {
      list <- mmMatrix[mmMatrix[,"construct"]==construct,"measurement"]
      
      for (item in list){
        adjustment <- adjustment + stats::sd(obsData[,item])*abs(as.numeric(outer_loadings[item,construct]))
        denom <- denom + abs(outer_loadings[item,construct])
      }
      adjustment <- adjustment/denom
      construct_scores[,construct] <- construct_scores[,construct]*adjustment
    }
  }
  return(construct_scores)
  
}

estimate_path_coef <- function(smMatrix, construct_scores,dependant, paths_matrix) {
  #Regression betas for the incoming paths
  #Iterate and regress the incoming paths
  for (i in 1:length(dependant))  {
    #Indentify the independant variables
    independant<-smMatrix[smMatrix[,"target"]==dependant[i],"source"]
    
    #Solve the system of equations
    paths_matrix[independant,dependant[i]] = solve(t(construct_scores[,independant]) %*% construct_scores[,independant], t(construct_scores[,independant]) %*% construct_scores[,dependant[i]])
  }
  return(paths_matrix)
}

# Calculate insample metrics ----
calc_insample <- function(obsData, construct_scores, smMatrix, dependant, construct_score_cors) {
  # matrix includes BIC
  # Remove BIC for now
  #insample <- matrix(, nrow=3, ncol=length(dependant), byrow =TRUE, dimnames = list(c("Rsq","AdjRsq","BIC"), dependant))
  
  # matrix excludes BIC
  insample <- matrix(, nrow=2, ncol=length(dependant), byrow =TRUE, dimnames = list(c("Rsq", "AdjRsq"), dependant))
  
  for (i in 1:length(dependant))  {
    #Indentify the independant variables
    independant <- smMatrix[smMatrix[, "target"]==dependant[i], "source"]
    
    #Calculate insample for endogenous
    #    construct_score_cors <- stats::cor(construct_scores)
    r_sq <- 1 - 1/solve(construct_score_cors[c(independant, dependant[i]), c(independant, dependant[i])])
    insample[1, i] <- r_sq[dependant[i], dependant[i]]
    insample[2, i] <- 1 - (1 - insample[1, i])*((nrow(obsData)-1)/(nrow(obsData)-length(independant) - 1))
    # Calculate the BIC for the endogenous
    # Remove BIC for now
    #insample[3, i] <- BIC_func(r_sq[dependant[i], dependant[i]], length(independant), nrow(obsData), construct_scores[, dependant[i]])
  }
  return(insample)
}

mean_replacement <- function(data) {
  for (i in 1:ncol(data)) {
    colmean <- mean(data[,i][!(is.na(data[,i]))])
    data[,i][is.na(data[,i])] <- colmean
  }
  return(data)
}


# HTMT as per Henseler, J., Ringle, C. M., & Sarstedt, M. (2014). A new criterion for assessing discriminant validity in
# variance-based structural equation modeling. Journal of the Academy of Marketing Science, 43(1), 115-135.
# https://doi.org/10.1007/s11747-014-0403-8
HTMT <- function(seminr_model) {
  if (is.null(seminr_model$hoc)) {
    constructs <- intersect(unique(seminr_model$smMatrix),unique(seminr_model$mmMatrix[,1 ]))
  } else {
    constructs <- intersect(unique(c(seminr_model$smMatrix, seminr_model$first_stage_model$smMatrix)),unique(seminr_model$mmMatrix[,1 ]))
  }
  
  HTMT <- matrix(, nrow=length(constructs), ncol=length(constructs),
                 dimnames = list(constructs,constructs))
  for (constructi in constructs[1:(length(constructs)-1)]) {
    for (constructj in constructs[(which(constructs == constructi)+1):length(constructs)]) {
      manifesti <- seminr_model$mmMatrix[seminr_model$mmMatrix[, 1] == constructi, "measurement"]
      manifestj <- seminr_model$mmMatrix[seminr_model$mmMatrix[, 1] == constructj, "measurement"]
      item_correlation_matrix <- abs(stats::cor(seminr_model$data[, manifesti],seminr_model$data[, manifestj]))
      HTHM <- mean(item_correlation_matrix)
      if(length(manifesti)>1 ) {
        cor_matrix <- abs(stats::cor(seminr_model$data[, manifesti], seminr_model$data[, manifesti]))
        diag(cor_matrix) <- 0
        MTHM <- (2/(length(manifesti)*(length(manifesti)-1)))*(sum(cor_matrix[!lower.tri(cor_matrix)]))
      } else {
        MTHM <- 1
      }
      if(length(manifestj)>1) {
        cor_matrix2 <- abs(stats::cor(seminr_model$data[, manifestj], seminr_model$data[, manifestj]))
        diag(cor_matrix2) <- 0
        MTHM <- sqrt(MTHM * (2/(length(manifestj)*(length(manifestj)-1)))*(sum(cor_matrix2[!lower.tri(cor_matrix2)])))
      } else {
        MTHM <- sqrt(1 * MTHM)
      }
      HTMT[constructi, constructj] <- HTHM / MTHM
    }
  }
  convert_to_table_output(HTMT)
}

# Function to return the total effects of a model
total_effects <- function(path_coef) {
  output <- path_coef
  paths <- path_coef
  while (sum(paths) > 0) {
    paths <- paths %*% path_coef
    output <- output + paths
  }
  return(output)
}

convert_to_table_output <- function(matrix) {
  class(matrix) <- append(class(matrix), "table_output")
  return(matrix)
}
