# FROM seminr


## original
cron_alpha <- function(cov_mat) {
  k <- nrow(cov_mat)
  cov_i <- sum(diag(cov_mat))
  alpha <- (k/(k-1))*(1 - (cov_i/sum(cov_mat)))
  return(alpha)
}

## original
cronbachs_alpha <- function(seminr_model) {
  constructs = colnames(plsModel$path_coef)
  alpha_vec <- c()
  for (i in constructs) {
    items <- seminr_model$mmMatrix[seminr_model$mmMatrix[,"construct"] == i,"measurement"]
    if (length(items) > 1) {
      cov_mat <- stats::cor(seminr_model$rawdata, seminr_model$rawdata)[items, items]
      alpha_vec[[i]] <- cron_alpha(cov_mat)
    } else {
      alpha_vec[[i]] <- 1
    }
  }
  return(unlist(alpha_vec))
}

# calculate_CR <- function(seminr_model) {
#   loadings <- seminr_model$outer_loadings
  
#   CRs <- sapply(1:ncol(loadings), function(i) {
#     lambda <- loadings[, i, drop = FALSE]
#     lambda_sq <- lambda^2
#     CR <- (sum(lambda))^2 / ((sum(lambda))^2 + sum(1 - lambda_sq))
#     CR
#   }, simplify = TRUE)
  
#   names(CRs) <- colnames(loadings)
#   return(CRs)
# }

## use comp_group to denote whether use the block of compositinal variable
calculate_CR <- function(seminr_model, comp_group = NULL) {
  
  if (is.null(comp_group)) {
    loadings <- rowSums(seminr_model$outer_loadings)
    # Calculate CR for entire set of manifest variables
    CRs <- sapply(1:ncol(seminr_model$outer_loadings), function(i) {
      lv <- colnames(seminr_model$outer_loadings)[i]
      mvs <- plsModel$mmMatrix[plsModel$mmMatrix[,"construct"] == lv, "measurement"]
      lambda <- loadings[mvs]
      lambda_sq <- lambda^2
      CR <- (sum(lambda))^2 / ((sum(lambda))^2 + sum(1 - lambda_sq))
      CR
    }, simplify = TRUE)
    names(CRs) <- colnames(seminr_model$outer_loadings)
    
  } else {
    loadings <- rowSums(seminr_model$outer_loadings)
    # Calculate CR for each subset defined by comp_group
    unique_groups <- unique(comp_group)
    CRs <- sapply(unique_groups, function(group) {
      lambda <- loadings[which(comp_group == group)]
      lambda_sq <- lambda^2
      CR <- (sum(lambda))^2 / ((sum(lambda))^2 + sum(1 - lambda_sq))
      CR
    }, simplify = TRUE)
    names(CRs) <- unique_groups
  }
  
  return(CRs)
}

# calculate_AVE <- function(seminr_model) {
#   loadings <- seminr_model$outer_loadings
  
#   AVEs <- sapply(1:ncol(loadings), function(i) {
#     lambda <- loadings[, i, drop = FALSE]
#     AVE <- sum(lambda^2) / length(lambda)
#     AVE
#   }, simplify = TRUE)
  
#   names(AVEs) <- colnames(loadings)
#   return(AVEs)
# }

## use comp_group to denote whether use the block of compositinal variable
calculate_AVE <- function(seminr_model, comp_group = NULL) {
  if (is.null(comp_group)) {
    loadings <- rowSums(seminr_model$outer_loadings)
    # Calculate AVE for entire set of manifest variables
    AVEs <- sapply(1:ncol(seminr_model$outer_loadings), function(i) {
      lv <- colnames(seminr_model$outer_loadings)[i]
      mvs <- plsModel$mmMatrix[plsModel$mmMatrix[,"construct"] == lv, "measurement"]
      lambda <- loadings[mvs]
      AVE <- sum(lambda^2) / length(lambda)
      AVE
    }, simplify = TRUE)
    names(AVEs) <- colnames(seminr_model$outer_loadings)

  } else {
    loadings <- rowSums(seminr_model$outer_loadings)
    # Calculate AVE for each subset defined by comp_group
    unique_groups <- unique(comp_group)
    AVEs <- sapply(unique_groups, function(group) {
      lambda <- loadings[which(comp_group == group)]
      AVE <- sum(lambda^2) / length(lambda)
      AVE
    }, simplify = TRUE)
    names(AVEs) <- unique_groups
  }
  
  return(AVEs)
}


# from seminr
## original
calculate_HTMT <- function(seminr_model) {
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

## original
calculate_R2 = function(seminr_model) {
    Rsquard = c()
    constructs = colnames(plsModel$path_coef)
    for (i in 1:length(constructs)) Rsquard = c(Rsquard, cor(seminr_model$construct_scores)[i,]%*%seminr_model$path_coef[,i])
    names(Rsquard) = constructs
    return(unlist(Rsquard))
}

## original
calculate_Q2 <- function(seminr_model, D = 7) {
  dataset = seminr_model$data
  construct_scores <- as.data.frame(seminr_model$construct_scores)
  loadings <- seminr_model$outer_loadings

  constructs <- names(construct_scores)
  
  Q2_values <- numeric(length(constructs))

  for (lv in constructs) {
    manifest_vars <- rownames(loadings)[loadings[,lv, drop=FALSE] != 0]
    
    SSE <- 0
    SSO <- 0
    
    for (var in manifest_vars) {
      actual_values <- dataset[,var]
      n <- length(actual_values)
      estimated_values <- construct_scores[[lv]] * loadings[var, lv]
      if (D < n) {
        # 使用遗漏距离 D 进行 Blindfolding
        omitted_indices <- seq(from = 1, to = n, by = D)
        omitted_values <- actual_values[omitted_indices]
        predicted_values <- estimated_values[omitted_indices]

        # 计算 SSE 和 SSO
        SSE <- SSE + sum((omitted_values - predicted_values)^2)
        SSO <- SSO + length(omitted_values)
      }
    }
    Q2 <- 1 - (SSE / SSO)
    Q2_values[which(constructs == lv)] <- Q2
  }
  
  names(Q2_values) <- constructs
  return(Q2_values)
}

## original
calculate_bootstrap <- function(seminr_model, nboot = 10000, cores = 2, seed) {
    # Bootstrapping for significance as per Hair, J. F., Hult, G. T. M., Ringle, C. M., and Sarstedt, M. (2017). A Primer on
    # Partial Least Squares Structural Equation Modeling (PLS-SEM), 2nd Ed., Sage: Thousand Oaks.
    message("Bootstrapping model using seminr...")
    
    # prepare parameters for cluster export (model parameters)
    d <- seminr_model$rawdata
    measurement_model <- seminr_model$measurement_model
    structural_model <- seminr_model$smMatrix
    inner_weights <- seminr_model$inner_weights
    # missing_value <- seminr_model$settings$missing_value
    missing_value <- NA
    # maxIt <- seminr_model$settings$maxIt
    # stopCriterion <- seminr_model$settings$stopCriterion
    # missing <- seminr_model$settings$missing
    missing <- mean_replacement
    
    
    # Initialize the cluster
    suppressWarnings(ifelse(is.null(cores), cl <- parallel::makeCluster(parallel::detectCores(), setup_strategy = "sequential"), cl <- parallel::makeCluster(cores, setup_strategy = "sequential")))
    
    # Function to generate random samples with replacement
    getRandomIndex <- function(d) {return(sample.int(nrow(d), replace = TRUE))}
    
    # Check for and create random seed if NULL
    if (is.null(seed)) {seed <- sample.int(100000, size = 1)}
    
    # Export variables and functions to cluster
    parallel::clusterExport(cl=cl, varlist=c("measurement_model",
                                            "structural_model",
                                            "inner_weights",
                                            "getRandomIndex",
                                            "d",
                                            "HTMT",
                                            "seed",
                                            "total_effects",
                                            "missing_value",
                                            "convert_to_table_output",
                                            "maxIt",
                                            "stopCriterion",
                                            "missing"), envir=environment())
    
    # Calculate the expected nrow of the bootmatrix
    length <- 3*nrow(seminr_model$path_coef)^2 + 2*nrow(seminr_model$outer_loadings)*ncol(seminr_model$outer_loadings)
    
    # Function to get PLS estimate results
    getEstimateResults <- function(i, d = d, length) {
        set.seed(seed + i)
        # plsc-bootstrap
        tryCatch({
        boot_model <- seminr::estimate_pls(data = d[getRandomIndex(d),],
                                            measurement_model,
                                            structural_model,
                                            inner_weights = inner_weights)
        boot_htmt <- HTMT(boot_model)
        boot_total <- total_effects(boot_model$path_coef)
        return(as.matrix(c(c(boot_model$path_coef),
                            c(boot_model$outer_loadings),
                            c(boot_model$outer_weights),
                            c(boot_htmt),
                            c(boot_total))))
        },
        error = function(cond) {
        message("Bootstrapping encountered an ERROR: ")
        message(cond)
        return(rep(NA, length))
        },
        warning = function(cond) {
        message("Bootstrapping encountered an ERROR: ")
        message(cond)
        return(rep(NA, length))
        }
        )
    }
    
    # Bootstrap the estimates
    # utils::capture.output(bootmatrix <- parallel::parSapply(cl, 1:nboot, getEstimateResults, d, length))
    bootmatrix <- sapply(1:nboot,getEstimateResults, d, length)
    
    # Clean the NAs and report the NAs
    bootmatrix <- bootmatrix[,!is.na(bootmatrix[1,])]
    fails <- nboot - ncol(bootmatrix)
    nboot <- nboot - fails
    if (fails > 0) {
        message(paste("Bootstrapping encountered a WARNING: ", fails, "bootstrap iterations failed to converge (possibly due to PLSc). \nThese failed iterations are excluded from the reported bootstrap statistics."))
    }
    
    
    # Collect means and sds for all estimates from bootmatrix
    means <- apply(bootmatrix,1,mean)
    sds <- apply(bootmatrix,1,stats::sd)
    
    # Create the matrix and array of bootstrapped paths ----
    
    # Collect the dimensions of the path_coef matrix
    path_cols <- ncol(seminr_model$path_coef)
    path_rows <- nrow(seminr_model$path_coef)
    
    # Identify start and end points for path coef data
    start <- 1
    end <- (path_cols*path_rows)
    
    # Create the array of bootstrap paths
    boot_paths <- array(bootmatrix[start:end,1:nboot], dim = c(path_rows, path_cols, nboot), dimnames = list(rownames(seminr_model$path_coef), colnames(seminr_model$path_coef),1:nboot))
    
    # Create the summary matrices of means and sds for path coefs
    paths_means <- matrix(means[start:end], nrow = path_rows, ncol = path_cols)
    paths_sds <- matrix(sds[start:end], nrow = path_rows, ncol = path_cols)
    
    # Create the bootstrapped paths matrix (take care to not be stranded with a single column/row vector)
    paths_descriptives <- cbind(seminr_model$path_coef, paths_means, paths_sds)
    
    # Clean the empty paths (take care to not be stranded with a single column/row vector)
    filled_cols <- apply(paths_descriptives != 0, 2, any, na.rm=TRUE)
    filled_rows <- apply(paths_descriptives != 0, 1, any, na.rm=TRUE)
    paths_descriptives <- subset(paths_descriptives, filled_rows, filled_cols)
    
    # Get the number of DVs
    if (length(unique(structural_model[,"target"])) == 1) {
        dependant <- unique(structural_model[,"target"])
    } else {
        dependant <- colnames(paths_descriptives[, 1:length(unique(structural_model[,"target"]))])
    }
    
    # Construct the vector of column names
    col_names <- c()
    # Clean the column names
    for (parameter in c("PLS Est.", "Boot Mean", "Boot SD")) {
        for (i in 1:length(dependant)) {
        col_names <- c(col_names, paste(dependant[i], parameter, sep = " "))
        }
    }
    
    # Assign column names
    colnames(paths_descriptives) <- col_names
    
    # Create the matrix and array of bootstrapped loadings ----
    
    # Collect the dimensions of the path_coef matrix
    mm_cols <- ncol(seminr_model$outer_loadings)
    mm_rows <- nrow(seminr_model$outer_loadings)
    
    # Identify start and end points for path coef data
    start <- end+1
    end <- start+(mm_cols*mm_rows)-1
    
    # create the array of bootstrapped loadings
    boot_loadings <- array(bootmatrix[start:end,1:nboot], dim = c(mm_rows, mm_cols, nboot), dimnames = list(rownames(seminr_model$outer_loadings), colnames(seminr_model$outer_loadings),1:nboot))
    
    # Create the summary matrices of means and sds for loadings
    loadings_means <- matrix(means[start:end], nrow = mm_rows, ncol = mm_cols)
    loadings_sds <- matrix(sds[start:end], nrow = mm_rows, ncol = mm_cols)
    
    # Subset paths matrix (take care to not be stranded with a single column/row vector)
    loadings_descriptives <- cbind(seminr_model$outer_loadings, loadings_means, loadings_sds)
    
    # Construct the vector of column names 2 for outer model
    col_names2 <- c()
    # Clean the column names
    for (parameter in c("PLS Est.", "Boot Mean", "Boot SD")) {
        for(i in colnames(seminr_model$outer_loadings)) {
        col_names2 <- c(col_names2, paste(i, parameter, sep = " "))
        }
    }
    
    # Assign column names to loadings matrix
    colnames(loadings_descriptives) <- col_names2
    
    # Identify start and end points for path weights data
    start <- end+1
    end <- start+(mm_cols*mm_rows)-1
    
    # create the array of bootstrapped weights
    boot_weights <- array(bootmatrix[start:end,1:nboot], dim = c(mm_rows, mm_cols, nboot), dimnames = list(rownames(seminr_model$outer_loadings), colnames(seminr_model$outer_loadings),1:nboot))
    
    # Create the summary matrices of means and sds
    weights_means <- matrix(means[start:end], nrow = mm_rows, ncol = mm_cols)
    weights_sds <- matrix(sds[start:end], nrow = mm_rows, ncol = mm_cols)
    
    # create the weights matrix (take care to not be stranded with a single column/row vector)
    weights_descriptives <- cbind(seminr_model$outer_weights, weights_means, weights_sds)
    
    # Assign column names to weights matrix
    colnames(weights_descriptives) <- col_names2
    
    # Collect HTMT matrix
    HTMT_matrix <- HTMT(seminr_model)
    
    # Identify start and end points for HTMT data
    start <- end+1
    end <- start+(ncol(HTMT_matrix)*nrow(HTMT_matrix))-1
    
    # Collect the array of bootstrapped HTMT
    boot_HTMT <- array(bootmatrix[start:end,1:nboot], dim = c(nrow(HTMT_matrix), ncol(HTMT_matrix), nboot), dimnames = list(rownames(HTMT_matrix), colnames(HTMT_matrix),1:nboot))
    
    # Collect the matrices of means and sds for HTMT
    htmt_means <- matrix(means[start:end], nrow = nrow(HTMT_matrix), ncol = ncol(HTMT_matrix))
    htmt_sds <- matrix(sds[start:end], nrow = nrow(HTMT_matrix), ncol = ncol(HTMT_matrix))
    
    # create HTMT matrix (take care to not be stranded with a single column/row vector)
    HTMT_descriptives <- cbind(HTMT_matrix, htmt_means, htmt_sds)
    
    # Construct the vector of column names 3 for HTMT
    col_names3 <- c()
    # Clean the column names
    for (parameter in c("PLS Est.", "Boot Mean", "Boot SD")) {
        for(i in colnames(HTMT_matrix)) {
        col_names3 <- c(col_names3, paste(i, parameter, sep = " "))
        }
    }
    # Get boot_HTMT column names
    colnames(HTMT_descriptives) <- col_names3
    
    # Subset total paths matrix (take care to not be stranded with a single column/row vector)
    total_matrix <- total_effects(seminr_model$path_coef)
    start <- end+1
    end <- start+(path_cols*path_rows)-1
    
    # Collect the array of bootstrapped total paths
    boot_total_paths <- array(bootmatrix[start:end,1:nboot], dim = c(path_rows, path_cols, nboot), dimnames = list(rownames(total_matrix), colnames(total_matrix),1:nboot))
    
    # Collect the matrices for means and sds
    total_means <- matrix(means[start:end], nrow = path_rows, ncol = path_cols)
    total_sds <- matrix(sds[start:end], nrow = path_rows, ncol = path_cols)
    
    
    # Subset paths matrix (take care to not be stranded with a single column/row vector)
    total_paths_descriptives <- cbind(total_matrix, total_means, total_sds)
    
    # Clean the empty paths (take care to not be stranded with a single column/row vector)
    filled_cols <- apply(total_paths_descriptives != 0, 2, any, na.rm=TRUE)
    filled_rows <- apply(total_paths_descriptives != 0, 1, any, na.rm=TRUE)
    total_paths_descriptives <- subset(total_paths_descriptives, filled_rows, filled_cols)
    
    # Assign column names
    colnames(total_paths_descriptives) <- col_names
    
    # Change the class of the descriptives objects
    class(paths_descriptives) <- append(class(paths_descriptives), "table_output")
    class(loadings_descriptives) <- append(class(loadings_descriptives), "table_output")
    class(weights_descriptives) <- append(class(weights_descriptives), "table_output")
    class(HTMT_descriptives) <- append(class(HTMT_descriptives), "table_output")
    class(total_paths_descriptives) <- append(class(total_paths_descriptives), "table_output")
    
    # Add the bootstrap matrix to the seminr_model object
    seminr_model$boot_paths <- boot_paths
    seminr_model$boot_loadings <- boot_loadings
    seminr_model$boot_weights <- boot_weights
    seminr_model$boot_HTMT <- boot_HTMT
    seminr_model$boot_total_paths <- boot_total_paths
    seminr_model$paths_descriptives <- paths_descriptives
    seminr_model$loadings_descriptives <- loadings_descriptives
    seminr_model$weights_descriptives <- weights_descriptives
    seminr_model$HTMT_descriptives <- HTMT_descriptives
    seminr_model$total_paths_descriptives <- total_paths_descriptives
    seminr_model$boots <- nboot
    seminr_model$seed <- seed
    class(seminr_model) <- c("boot_seminr_model", "seminr_model")
    message("SEMinR Model successfully bootstrapped")
    parallel::stopCluster(cl)
    bootmodel <- seminr_model

    return(summary(bootmodel)[[2]])

}