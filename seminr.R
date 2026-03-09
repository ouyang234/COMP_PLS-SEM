# ==============================================================================
# seminr.R - Custom PLS-SEM modeling functions
# ==============================================================================

# --- Measurement Model Definition ---

constructs <- function(...) {
  return_list <- list(...)
  names(return_list) <- lapply(return_list, function(x) class(x)[[3]])
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
    composite_type <- "A"
  } else if (identical(weights, regression_weights)) {
    composite_type <- "B"
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

# --- Weighting Schemes ---

mode_A <- function(mmMatrix, i, normData, construct_scores) {
  return(stats::cov(
    normData[, mmMatrix[mmMatrix[, "construct"] == i, "measurement"]],
    construct_scores[, i]
  ))
}
correlation_weights <- mode_A

mode_B <- function(mmMatrix, i, normData, construct_scores) {
  items <- mmMatrix[mmMatrix[, "construct"] == i, "measurement"]
  return(
    solve(stats::cor(normData[, items])) %*%
      stats::cor(normData[, items], construct_scores[, i])
  )
}
regression_weights <- mode_B

# --- Structural Model Definition ---

relationships <- function(...) {
  smMatrix <- matrix(c(...), ncol = 2, byrow = TRUE,
                     dimnames = list(NULL, c("source", "target")))
  return(smMatrix)
}

paths <- function(from, to) {
  return(as.vector(t(as.matrix(expand.grid(from, to)))))
}

# --- Helper Functions ---

mmMatrix_per_construct <- function(construct, mmMatrix) {
  constructmatrix <- mmMatrix[mmMatrix[, "construct"] == construct,
                              c("construct", "measurement", "type")]
  if (class(constructmatrix)[1] != "matrix") {
    constructmatrix <- t(as.matrix(constructmatrix))
  }
  return(constructmatrix)
}

compute_construct_rhoA <- function(weights, mmMatrix, construct, obsData) {
  w <- as.matrix(weights[mmMatrix[mmMatrix[, "construct"] == construct, "measurement"], construct])
  indicators <- scale(obsData[, mmMatrix[mmMatrix[, "construct"] == construct, "measurement"]], TRUE, TRUE)
  S <- stats::cov(indicators, indicators)
  diag(S) <- 0

  AAnondiag <- w %*% t(w)
  diag(AAnondiag) <- 0

  return((t(w) %*% w)^2 * ((t(w) %*% S %*% w) / (t(w) %*% AAnondiag %*% w)))
}

all_reflective <- function(mmMatrix, constructs) {
  unique(mmMatrix[mmMatrix[, "type"] == "C", "construct"])
}

interaction_term <- function(iv, moderator, method = product_indicator, weights = mode_A) {
  intxn <- method(iv, moderator, weights)
  class(intxn) <- class(method())
  return(intxn)
}

construct_names <- function(smMatrix) {
  unique(c(smMatrix[, 1], smMatrix[, 2]))
}

all_endogenous <- function(smMatrix) {
  unique(smMatrix[, "target"])
}

get_measure_mode <- function(construct, mmMatrix) {
  type <- mmMatrix[mmMatrix[, "construct"] == construct, "type"][1]
  if (type %in% c("A", "C", "HOCA")) {
    return(mode_A)
  } else if (type %in% c("B", "HOCB")) {
    return(mode_B)
  } else if (type == "UNIT") {
    return(unit_weights)
  }
}

# --- Standardization & Data Preparation ---

standardize_safely <- function(x) {
  res <- scale(x, TRUE, TRUE)
  if (any(is.nan(res))) stop("zero variance items cannot be scaled")
  res
}

standardize_outer_weights <- function(normData, mmVariables, outer_weights) {
  std_devs <- attr(scale(normData[, mmVariables] %*% outer_weights, center = FALSE), "scaled:scale")
  return(t(t(outer_weights) / std_devs))
}

calculate_loadings <- function(weights_matrix, construct_scores, normData) {
  return(as.matrix(stats::cov(normData, construct_scores) * weights_matrix))
}

mean_replacement <- function(data) {
  for (i in 1:ncol(data)) {
    colmean <- mean(data[, i], na.rm = TRUE)
    data[is.na(data[, i]), i] <- colmean
  }
  return(data)
}

# --- Inner Weighting Schemes ---

path_factorial <- function(smMatrix, construct_scores, dependant, paths_matrix) {
  inner_paths <- stats::cor(construct_scores, construct_scores) * (paths_matrix + t(paths_matrix))
  return(inner_paths)
}

path_weighting <- function(smMatrix, construct_scores, dependant, paths_matrix) {
  inner_paths <- stats::cor(construct_scores, construct_scores) * t(paths_matrix)

  for (i in seq_along(dependant)) {
    independant <- smMatrix[smMatrix[, "target"] == dependant[i], "source"]
    inner_paths[independant, dependant[i]] <- solve(
      t(construct_scores[, independant]) %*% construct_scores[, independant],
      t(construct_scores[, independant]) %*% construct_scores[, dependant[i]]
    )
  }
  return(inner_paths)
}

# --- Interaction Term Adjustment ---
# Henseler & Chin (2010), Structural Equation Modeling, 17(1), 82-109.
adjust_interaction <- function(constructs, mmMatrix, outer_loadings, construct_scores, obsData) {
  for (construct in constructs) {
    if (grepl("\\*", construct)) {
      items <- mmMatrix[mmMatrix[, "construct"] == construct, "measurement"]
      adjustment <- 0
      denom <- 0
      for (item in items) {
        adjustment <- adjustment + stats::sd(obsData[, item]) * abs(as.numeric(outer_loadings[item, construct]))
        denom <- denom + abs(outer_loadings[item, construct])
      }
      construct_scores[, construct] <- construct_scores[, construct] * (adjustment / denom)
    }
  }
  return(construct_scores)
}

# --- Path Coefficient Estimation ---

estimate_path_coef <- function(smMatrix, construct_scores, dependant, paths_matrix) {
  for (i in seq_along(dependant)) {
    independant <- smMatrix[smMatrix[, "target"] == dependant[i], "source"]
    paths_matrix[independant, dependant[i]] <- solve(
      t(construct_scores[, independant]) %*% construct_scores[, independant],
      t(construct_scores[, independant]) %*% construct_scores[, dependant[i]]
    )
  }
  return(paths_matrix)
}

# --- In-Sample Metrics ---

calc_insample <- function(obsData, construct_scores, smMatrix, dependant, construct_score_cors) {
  insample <- matrix(, nrow = 2, ncol = length(dependant), byrow = TRUE,
                     dimnames = list(c("Rsq", "AdjRsq"), dependant))

  for (i in seq_along(dependant)) {
    independant <- smMatrix[smMatrix[, "target"] == dependant[i], "source"]
    r_sq <- 1 - 1 / solve(construct_score_cors[c(independant, dependant[i]),
                                                 c(independant, dependant[i])])
    insample[1, i] <- r_sq[dependant[i], dependant[i]]
    insample[2, i] <- 1 - (1 - insample[1, i]) *
      ((nrow(obsData) - 1) / (nrow(obsData) - length(independant) - 1))
  }
  return(insample)
}

# --- HTMT ---
# Henseler, Ringle & Sarstedt (2014), JAMS, 43(1), 115-135.
HTMT <- function(seminr_model) {
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

      item_correlation_matrix <- abs(stats::cor(seminr_model$data[, manifesti],
                                                seminr_model$data[, manifestj]))
      HTHM <- mean(item_correlation_matrix)

      if (length(manifesti) > 1) {
        cor_matrix <- abs(stats::cor(seminr_model$data[, manifesti], seminr_model$data[, manifesti]))
        diag(cor_matrix) <- 0
        MTHM <- (2 / (length(manifesti) * (length(manifesti) - 1))) *
          sum(cor_matrix[!lower.tri(cor_matrix)])
      } else {
        MTHM <- 1
      }

      if (length(manifestj) > 1) {
        cor_matrix2 <- abs(stats::cor(seminr_model$data[, manifestj], seminr_model$data[, manifestj]))
        diag(cor_matrix2) <- 0
        MTHM <- sqrt(MTHM * (2 / (length(manifestj) * (length(manifestj) - 1))) *
                       sum(cor_matrix2[!lower.tri(cor_matrix2)]))
      } else {
        MTHM <- sqrt(MTHM)
      }

      HTMT[constructi, constructj] <- HTHM / MTHM
    }
  }
  convert_to_table_output(HTMT)
}

# --- Total Effects ---

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
