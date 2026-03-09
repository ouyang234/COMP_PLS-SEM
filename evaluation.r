# ==============================================================================
# evaluation.r - Model evaluation functions (CR, AVE, HTMT, R2, Q2, Bootstrap)
# ==============================================================================

# --- Cronbach's Alpha ---

cron_alpha <- function(cov_mat) {
  k <- nrow(cov_mat)
  cov_i <- sum(diag(cov_mat))
  (k / (k - 1)) * (1 - cov_i / sum(cov_mat))
}

cronbachs_alpha <- function(seminr_model) {
  constructs <- colnames(seminr_model$path_coef)
  alpha_vec <- c()
  for (i in constructs) {
    items <- seminr_model$mmMatrix[seminr_model$mmMatrix[, "construct"] == i, "measurement"]
    if (length(items) > 1) {
      cov_mat <- stats::cor(seminr_model$rawdata, seminr_model$rawdata)[items, items]
      alpha_vec[[i]] <- cron_alpha(cov_mat)
    } else {
      alpha_vec[[i]] <- 1
    }
  }
  return(unlist(alpha_vec))
}

# --- Composite Reliability ---
# Use comp_group to calculate CR at the compositional block level

calculate_CR <- function(seminr_model, comp_group = NULL) {
  loadings <- rowSums(seminr_model$outer_loadings)

  if (is.null(comp_group)) {
    CRs <- sapply(colnames(seminr_model$outer_loadings), function(lv) {
      mvs <- seminr_model$mmMatrix[seminr_model$mmMatrix[, "construct"] == lv, "measurement"]
      lambda <- loadings[mvs]
      sum(lambda)^2 / (sum(lambda)^2 + sum(1 - lambda^2))
    })
  } else {
    unique_groups <- unique(comp_group)
    CRs <- sapply(unique_groups, function(group) {
      lambda <- loadings[which(comp_group == group)]
      sum(lambda)^2 / (sum(lambda)^2 + sum(1 - lambda^2))
    })
    names(CRs) <- unique_groups
  }
  return(CRs)
}

# --- Average Variance Extracted ---
# Use comp_group to calculate AVE at the compositional block level

calculate_AVE <- function(seminr_model, comp_group = NULL) {
  loadings <- rowSums(seminr_model$outer_loadings)

  if (is.null(comp_group)) {
    AVEs <- sapply(colnames(seminr_model$outer_loadings), function(lv) {
      mvs <- seminr_model$mmMatrix[seminr_model$mmMatrix[, "construct"] == lv, "measurement"]
      lambda <- loadings[mvs]
      sum(lambda^2) / length(lambda)
    })
  } else {
    unique_groups <- unique(comp_group)
    AVEs <- sapply(unique_groups, function(group) {
      lambda <- loadings[which(comp_group == group)]
      sum(lambda^2) / length(lambda)
    })
    names(AVEs) <- unique_groups
  }
  return(AVEs)
}

# --- HTMT ---
# Henseler, Ringle & Sarstedt (2014), JAMS, 43(1), 115-135.

calculate_HTMT <- function(seminr_model) {
  if (is.null(seminr_model$hoc)) {
    constructs <- intersect(unique(seminr_model$smMatrix), unique(seminr_model$mmMatrix[, 1]))
  } else {
    constructs <- intersect(
      unique(c(seminr_model$smMatrix, seminr_model$first_stage_model$smMatrix)),
      unique(seminr_model$mmMatrix[, 1])
    )
  }

  HTMT <- matrix(, nrow = length(constructs), ncol = length(constructs),
                 dimnames = list(constructs, constructs))

  for (constructi in constructs[1:(length(constructs) - 1)]) {
    for (constructj in constructs[(which(constructs == constructi) + 1):length(constructs)]) {
      manifesti <- seminr_model$mmMatrix[seminr_model$mmMatrix[, 1] == constructi, "measurement"]
      manifestj <- seminr_model$mmMatrix[seminr_model$mmMatrix[, 1] == constructj, "measurement"]

      HTHM <- mean(abs(stats::cor(seminr_model$data[, manifesti], seminr_model$data[, manifestj])))

      if (length(manifesti) > 1) {
        cor_i <- abs(stats::cor(seminr_model$data[, manifesti], seminr_model$data[, manifesti]))
        diag(cor_i) <- 0
        MTHM <- (2 / (length(manifesti) * (length(manifesti) - 1))) * sum(cor_i[!lower.tri(cor_i)])
      } else {
        MTHM <- 1
      }

      if (length(manifestj) > 1) {
        cor_j <- abs(stats::cor(seminr_model$data[, manifestj], seminr_model$data[, manifestj]))
        diag(cor_j) <- 0
        MTHM <- sqrt(MTHM * (2 / (length(manifestj) * (length(manifestj) - 1))) *
                       sum(cor_j[!lower.tri(cor_j)]))
      } else {
        MTHM <- sqrt(MTHM)
      }

      HTMT[constructi, constructj] <- HTHM / MTHM
    }
  }
  convert_to_table_output(HTMT)
}

# --- R-squared ---

calculate_R2 <- function(seminr_model) {
  constructs <- colnames(seminr_model$path_coef)
  cor_scores <- cor(seminr_model$construct_scores)
  Rsquared <- sapply(seq_along(constructs), function(i) {
    cor_scores[i, ] %*% seminr_model$path_coef[, i]
  })
  names(Rsquared) <- constructs
  return(Rsquared)
}

# --- Q-squared (Blindfolding) ---

calculate_Q2 <- function(seminr_model, D = 7) {
  dataset <- seminr_model$data
  construct_scores <- as.data.frame(seminr_model$construct_scores)
  loadings <- seminr_model$outer_loadings
  constructs <- names(construct_scores)

  Q2_values <- numeric(length(constructs))

  for (lv in constructs) {
    manifest_vars <- rownames(loadings)[loadings[, lv, drop = FALSE] != 0]
    SSE <- 0
    SSO <- 0

    for (var in manifest_vars) {
      actual_values <- dataset[, var]
      n <- length(actual_values)
      estimated_values <- construct_scores[[lv]] * loadings[var, lv]

      if (D < n) {
        omitted_indices <- seq(from = 1, to = n, by = D)
        SSE <- SSE + sum((actual_values[omitted_indices] - estimated_values[omitted_indices])^2)
        SSO <- SSO + length(omitted_indices)
      }
    }
    Q2_values[which(constructs == lv)] <- 1 - SSE / SSO
  }

  names(Q2_values) <- constructs
  return(Q2_values)
}

# --- Bootstrap ---
# Hair, Hult, Ringle & Sarstedt (2017), A Primer on PLS-SEM, 2nd Ed., Sage.

calculate_bootstrap <- function(seminr_model, nboot = 10000, cores = 2, seed = NULL) {
  message("Bootstrapping model using seminr...")

  d <- seminr_model$rawdata
  measurement_model <- seminr_model$measurement_model
  structural_model  <- seminr_model$smMatrix
  inner_weights     <- seminr_model$inner_weights
  missing_value     <- NA
  missing           <- mean_replacement

  suppressWarnings(
    cl <- parallel::makeCluster(
      ifelse(is.null(cores), parallel::detectCores(), cores),
      setup_strategy = "sequential"
    )
  )

  if (is.null(seed)) seed <- sample.int(100000, size = 1)

  parallel::clusterExport(
    cl = cl,
    varlist = c("measurement_model", "structural_model", "inner_weights",
                "d", "HTMT", "seed", "total_effects", "missing_value",
                "convert_to_table_output", "maxIt", "stopCriterion", "missing"),
    envir = environment()
  )

  path_rows <- nrow(seminr_model$path_coef)
  path_cols <- ncol(seminr_model$path_coef)
  mm_rows   <- nrow(seminr_model$outer_loadings)
  mm_cols   <- ncol(seminr_model$outer_loadings)

  expected_length <- 3 * path_rows^2 + 2 * mm_rows * mm_cols

  getEstimateResults <- function(i, d, length) {
    set.seed(seed + i)
    tryCatch({
      boot_idx <- sample.int(nrow(d), replace = TRUE)
      boot_model <- seminr::estimate_pls(
        data = d[boot_idx, ],
        measurement_model, structural_model,
        inner_weights = inner_weights
      )
      boot_htmt  <- HTMT(boot_model)
      boot_total <- total_effects(boot_model$path_coef)
      as.matrix(c(c(boot_model$path_coef), c(boot_model$outer_loadings),
                   c(boot_model$outer_weights), c(boot_htmt), c(boot_total)))
    },
    error = function(cond) {
      message("Bootstrapping encountered an ERROR: ", cond)
      rep(NA, length)
    },
    warning = function(cond) {
      message("Bootstrapping encountered a WARNING: ", cond)
      rep(NA, length)
    })
  }

  bootmatrix <- sapply(1:nboot, getEstimateResults, d, expected_length)

  # Remove failed iterations
  bootmatrix <- bootmatrix[, !is.na(bootmatrix[1, ])]
  fails <- nboot - ncol(bootmatrix)
  nboot <- nboot - fails
  if (fails > 0) {
    message(fails, " bootstrap iterations failed to converge and are excluded.")
  }

  means <- apply(bootmatrix, 1, mean)
  sds   <- apply(bootmatrix, 1, stats::sd)

  # --- Extract bootstrapped paths ---
  start <- 1
  end   <- path_cols * path_rows

  boot_paths <- array(bootmatrix[start:end, 1:nboot],
                      dim = c(path_rows, path_cols, nboot),
                      dimnames = list(rownames(seminr_model$path_coef),
                                      colnames(seminr_model$path_coef), 1:nboot))

  paths_descriptives <- cbind(
    seminr_model$path_coef,
    matrix(means[start:end], nrow = path_rows, ncol = path_cols),
    matrix(sds[start:end],   nrow = path_rows, ncol = path_cols)
  )

  filled_cols <- apply(paths_descriptives != 0, 2, any, na.rm = TRUE)
  filled_rows <- apply(paths_descriptives != 0, 1, any, na.rm = TRUE)
  paths_descriptives <- subset(paths_descriptives, filled_rows, filled_cols)

  dependant <- if (length(unique(structural_model[, "target"])) == 1) {
    unique(structural_model[, "target"])
  } else {
    colnames(paths_descriptives[, 1:length(unique(structural_model[, "target"]))])
  }

  col_names <- as.vector(outer(dependant, c("PLS Est.", "Boot Mean", "Boot SD"), paste))
  colnames(paths_descriptives) <- col_names

  # --- Extract bootstrapped loadings ---
  start <- end + 1
  end   <- start + mm_cols * mm_rows - 1

  boot_loadings <- array(bootmatrix[start:end, 1:nboot],
                         dim = c(mm_rows, mm_cols, nboot),
                         dimnames = list(rownames(seminr_model$outer_loadings),
                                         colnames(seminr_model$outer_loadings), 1:nboot))

  loadings_descriptives <- cbind(
    seminr_model$outer_loadings,
    matrix(means[start:end], nrow = mm_rows, ncol = mm_cols),
    matrix(sds[start:end],   nrow = mm_rows, ncol = mm_cols)
  )

  col_names2 <- as.vector(outer(colnames(seminr_model$outer_loadings),
                                c("PLS Est.", "Boot Mean", "Boot SD"), paste))
  colnames(loadings_descriptives) <- col_names2

  # --- Extract bootstrapped weights ---
  start <- end + 1
  end   <- start + mm_cols * mm_rows - 1

  boot_weights <- array(bootmatrix[start:end, 1:nboot],
                        dim = c(mm_rows, mm_cols, nboot),
                        dimnames = list(rownames(seminr_model$outer_loadings),
                                        colnames(seminr_model$outer_loadings), 1:nboot))

  weights_descriptives <- cbind(
    seminr_model$outer_weights,
    matrix(means[start:end], nrow = mm_rows, ncol = mm_cols),
    matrix(sds[start:end],   nrow = mm_rows, ncol = mm_cols)
  )
  colnames(weights_descriptives) <- col_names2

  # --- Extract bootstrapped HTMT ---
  HTMT_matrix <- HTMT(seminr_model)
  start <- end + 1
  end   <- start + ncol(HTMT_matrix) * nrow(HTMT_matrix) - 1

  boot_HTMT <- array(bootmatrix[start:end, 1:nboot],
                     dim = c(nrow(HTMT_matrix), ncol(HTMT_matrix), nboot),
                     dimnames = list(rownames(HTMT_matrix), colnames(HTMT_matrix), 1:nboot))

  HTMT_descriptives <- cbind(
    HTMT_matrix,
    matrix(means[start:end], nrow = nrow(HTMT_matrix), ncol = ncol(HTMT_matrix)),
    matrix(sds[start:end],   nrow = nrow(HTMT_matrix), ncol = ncol(HTMT_matrix))
  )

  col_names3 <- as.vector(outer(colnames(HTMT_matrix), c("PLS Est.", "Boot Mean", "Boot SD"), paste))
  colnames(HTMT_descriptives) <- col_names3

  # --- Extract bootstrapped total effects ---
  total_matrix <- total_effects(seminr_model$path_coef)
  start <- end + 1
  end   <- start + path_cols * path_rows - 1

  boot_total_paths <- array(bootmatrix[start:end, 1:nboot],
                            dim = c(path_rows, path_cols, nboot),
                            dimnames = list(rownames(total_matrix), colnames(total_matrix), 1:nboot))

  total_paths_descriptives <- cbind(
    total_matrix,
    matrix(means[start:end], nrow = path_rows, ncol = path_cols),
    matrix(sds[start:end],   nrow = path_rows, ncol = path_cols)
  )

  filled_cols <- apply(total_paths_descriptives != 0, 2, any, na.rm = TRUE)
  filled_rows <- apply(total_paths_descriptives != 0, 1, any, na.rm = TRUE)
  total_paths_descriptives <- subset(total_paths_descriptives, filled_rows, filled_cols)
  colnames(total_paths_descriptives) <- col_names

  # Set table_output class on all descriptive matrices
  for (obj_name in c("paths_descriptives", "loadings_descriptives",
                     "weights_descriptives", "HTMT_descriptives",
                     "total_paths_descriptives")) {
    obj <- get(obj_name)
    class(obj) <- append(class(obj), "table_output")
    assign(obj_name, obj)
  }

  # Assemble output model
  seminr_model$boot_paths              <- boot_paths
  seminr_model$boot_loadings           <- boot_loadings
  seminr_model$boot_weights            <- boot_weights
  seminr_model$boot_HTMT               <- boot_HTMT
  seminr_model$boot_total_paths        <- boot_total_paths
  seminr_model$paths_descriptives      <- paths_descriptives
  seminr_model$loadings_descriptives   <- loadings_descriptives
  seminr_model$weights_descriptives    <- weights_descriptives
  seminr_model$HTMT_descriptives       <- HTMT_descriptives
  seminr_model$total_paths_descriptives <- total_paths_descriptives
  seminr_model$boots <- nboot
  seminr_model$seed  <- seed

  class(seminr_model) <- c("boot_seminr_model", "seminr_model")
  message("SEMinR Model successfully bootstrapped")
  parallel::stopCluster(cl)

  return(summary(seminr_model)[[2]])
}
