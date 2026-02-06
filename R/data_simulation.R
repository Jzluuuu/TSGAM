################################################################################
#' Gaussian distribution Data Simulation
#'
#' This function simulates data with margin- and nonlinear-dependent age effects
#' from Gaussian distribution for single-omics data.
#'
#' @importFrom stats rnorm model.matrix lm
#'
#' @param N Integer, the number of samples.
#' @param mean.pat Numeric, the mean of different pattern effect coefficients
#' for all related variables.
#' @param sd.pat Numeric, the variability (standard deviation) of different
#' pattern effect coefficients for all related variables.
#' @param mean.diff Numeric, the mean of different expressed effect coefficients
#' for all related variables
#' @param sd.diff Numeric, the variability (standard deviation) of different
#' expressed effect coefficients for all related variables.
#' @param mean.bg Numeric, the mean of background noise.
#' @param sd.bg Numeric, the variability (standard deviation) of
#' background noise.
#' @param p_total Integer, the number of variables.
#' @param p_pattern_relevant Integer, the number of variables affected by
#' different pattern effects.
#' @param p_diff_relevant Integer, the number of variables affected by
#' different expressed effects.
#' @param percentage_overlap_variables Numeric, in the range of \eqn{0}
#' to \eqn{1}, the percentage of variables with both confounder and treatment
#' effects given fixed number of variables with each effect.
#' @param seeds Integer. Default value is randomly selected from
#' \eqn{0} to \eqn{999}.
#'
#'
#' @return \code{simData} returns a list that contains the following components:
#'
#' \item{data}{The simulated data.}
#' \item{Age}{The simulated outcome variable indicating the Age information.}
#' \item{Group}{The simulated outcome variable indicating the Age group
#' information.}
#' \item{true.pat}{The variables assigned with the different pattern effect.}
#' \item{true.diff}{The variables assigned with the different expressed effect.}
#'
#' @author Jianzhong Lu
#'
#'
#' @examples
#' sim <- simData(N = 40,
#'                mean.pat = 2,
#'                sd.pat = 1,
#'                mean.diff = 2,
#'                sd.diff = 1,
#'                mean.bg = 0,
#'                sd.bg = 0.5,
#'                p_total = 2000,
#'                p_pattern_relevant = 40,
#'                p_diff_relevant = 40,
#'                percentage_overlap_variables = 0.5,
#'                seeds = 123)
#'
#'
#' @export
simData <- function(N = 40,
                    mean.pat = 2,
                    sd.pat = 1,
                    mean.diff = 0,
                    sd.diff = 1,
                    mean.bg = 0,
                    sd.bg = 0.5,
                    p_total = 500,
                    p_pattern_relevant = 60,
                    p_diff_relevant = 100,
                    percentage_overlap_variables = 0.5,
                    seeds = 123){

  set.seed(seeds)

  ### age info
  Young_age <- round(runif(N/2, min = 20, max = 35))
  Old_age   <- round(runif(N/2, min = 50, max = 65))
  Age <- c(Young_age, Old_age)
  Group <- rep(c("Young","Old"), each = N/2)

  ### design matrix (pattern): shape diff, remove mean diff
  YoungSeg <- as.integer(Group == "Young")
  OldSeg   <- as.integer(Group == "Old")

  uY <- (Age - min(Young_age)) / (max(Young_age) - min(Young_age) + 1e-8)
  uO <- (Age - min(Old_age))   / (max(Old_age)   - min(Old_age)   + 1e-8)
  uY[OldSeg == 1] <- 0
  uO[YoungSeg == 1] <- 0

  X_raw <- cbind(
    YoungSeg * uY,                        # Young linear trend
    OldSeg   * uO,                        # Old linear trend
    YoungSeg * cos(2*pi*uY),              # Young nonlinear
    OldSeg   * ((uO - 0.5)^2 - mean((uO[OldSeg==1]-0.5)^2))  # Old nonlinear
  )


  # residualize against Group to eliminate mean shift
  Z <- model.matrix(~ Group)  # (Intercept + GroupOld)
  X <- apply(X_raw, 2, function(v) as.numeric(residuals(lm(v ~ Z - 1))))

  ### pattern effects coefficients
  index.pat <- 1:p_pattern_relevant
  K_pat <- ncol(X)

  beta_pat <- matrix(0, nrow = p_total, ncol = K_pat)
  beta_pat[index.pat,] <- matrix(
    rnorm(K_pat * length(index.pat), mean = mean.pat, sd = sd.pat),
    ncol = K_pat
  )

  ### diff effects coefficients
  min.num <- min(p_pattern_relevant, p_diff_relevant)
  set.seed(999)
  index.inter <- sample(index.pat, min.num * percentage_overlap_variables)
  index.diff.only <- sample((1:p_total)[-index.pat],
                            (p_diff_relevant - min.num * percentage_overlap_variables))
  index.diff <- c(index.inter, index.diff.only)

  beta_diff <- matrix(0, nrow = p_total, ncol = 1)
  beta_diff[index.diff, 1] <- rnorm(p_diff_relevant, mean = mean.diff, sd = sd.diff)

  # build mu matrix
  G <- matrix(as.integer(Group == "Old"), ncol = 1)  # pure group effect
  mu.simu <- G %*% t(beta_diff) + X %*% t(beta_pat)

  # add noise
  noise.matrix <- matrix(rnorm(N*p_total, mean = mean.bg, sd = sd.bg), nrow = N)
  mu.data <- mu.simu + noise.matrix

  colnames(mu.data) <- paste0("F",1:p_total)

  return(list(data = mu.data,
              Age = Age,
              Group = Group,
              true.pat = index.pat,
              true.diff = index.diff))
}


################################################################################
#' Multi-Omics Data based on Gaussian distribution Simulation
#'
#' This function simulates data with margin- and nonlinear-dependent age effects
#' from Gaussian distribution for multi-omics data.
#'
#' @importFrom stats rnorm runif model.matrix
#'
#' @param N Integer, the number of samples.
#' @param mean.pat Numeric, the mean of different pattern effect coefficients
#' for all related variables.
#' @param sd.pat Numeric, the variability (standard deviation) of different
#' pattern effect coefficients for all related variables.
#' @param mean.diff Numeric, the mean of different expressed effect coefficients
#' for all related variables
#' @param sd.diff Numeric, the variability (standard deviation) of different
#' expressed effect coefficients for all related variables.
#' @param mean.bg Numeric, the mean of background noise.
#' @param sd.bg Numeric, the variability (standard deviation) of
#' background noise.
#' @param blocks List, including integers, the numbers of different variables.
#' @param p_pattern_relevant List, including integers, the numbers of different
#' variables affected by different pattern effects.
#' @param p_diff_relevant List, including integers, the numbers of different
#' variables affected by different expressed effects.
#' @param n_shared_factors Integer. Number of shared factors across omics
#' blocks. Larger values introduce stronger shared structure between
#' different omics layers. Set to \eqn{0} to disable shared effects.
#' @param shared_factor_sd Numeric. Standard deviation of the shared factors.
#' Controls the overall strength of the shared signal across omics blocks.
#' @param loading_sd Numeric. Standard deviation of the loadings linking shared
#' factors to individual features. Larger values increase feature-level
#' variation induced by shared factors.
#' @param percentage_overlap_variables Numeric, in the range of \eqn{0}
#' to \eqn{1}, the percentage of variables with both confounder and treatment
#' effects given fixed number of variables with each effect.
#' @param seeds Integer. Default value is randomly selected from
#' \eqn{0} to \eqn{999}.
#'
#'
#' @return \code{simData} returns a list that contains the following components:
#'
#' \item{data}{The simulated data.}
#' \item{Age}{The simulated outcome variable indicating the Age information.}
#' \item{Group}{The simulated outcome variable indicating the Age group
#' information.}
#' \item{true.pat}{The variables assigned with the different pattern effect.}
#' \item{true.diff}{The variables assigned with the different expressed effect.}
#'
#' @author Jianzhong Lu
#'
#'
#' @examples
#'sim <- simulate_multiomics(N = 40,
#'                           blocks = list(transcriptome = 2000,
#'                                         metabolome = 1000),
#'                           p_pattern_relevant = c(transcriptome = 40, metabolome = 20),
#'                           p_diff_relevant    = c(transcriptome = 40, metabolome = 20),
#'                           n_shared_factors = 3,
#'                           shared_factor_sd = 1,
#'                           loading_sd = 1,
#'                           mean.pat = 2, sd.pat = 1,
#'                           mean.diff = 2, sd.diff = 1,
#'                           mean.bg = 0, sd.bg = 0.5,
#'                           percentage_overlap_variables = 0.5,
#'                           seeds = 123)
#'
#'
#' @export
simulate_multiomics <- function(N = 40,
                                blocks = list(transcriptome = 500,
                                              proteome = 200),
                                p_pattern_relevant = c(transcriptome = 20, proteome = 8),
                                p_diff_relevant    = c(transcriptome = 40, proteome = 15),
                                n_shared_factors = 3,
                                shared_factor_sd = 1,
                                loading_sd = 1,
                                mean.pat = 2, sd.pat = 1,
                                mean.diff = 0, sd.diff = 1,
                                mean.bg = 0, sd.bg = 0.5,
                                percentage_overlap_variables = 0.5,
                                seeds = 123
) {
  set.seed(seeds)

  if (!is.list(blocks) || is.null(names(blocks)) || any(names(blocks) == "")) {
    stop("`blocks` must be a named list, e.g., list(transcriptome=500, proteome=200).")
  }
  if (N %% 2 != 0) stop("N must be even (Young/Old each N/2).")

  # -----------------------------
  # Sample-level phenotype
  # -----------------------------
  Young_age <- round(runif(N/2, min = 20, max = 35))
  Old_age   <- round(runif(N/2, min = 50, max = 65))
  Age <- c(Young_age, Old_age)
  Group <- rep(c("Young", "Old"), each = N/2)

  YoungSeg <- as.integer(Group == "Young")
  OldSeg   <- as.integer(Group == "Old")

  # normalized age within group
  uY <- (Age - min(Young_age)) / (max(Young_age) - min(Young_age) + 1e-8)
  uO <- (Age - min(Old_age))   / (max(Old_age)   - min(Old_age)   + 1e-8)
  uY[OldSeg == 1] <- 0
  uO[YoungSeg == 1] <- 0

  # pattern design matrix
  X_raw <- cbind(
    YoungSeg * uY,
    OldSeg   * uO,
    YoungSeg * cos(2*pi*uY),
    OldSeg   * ((uO - 0.5)^2 - mean((uO[OldSeg == 1] - 0.5)^2))
  )

  # remove group mean shift
  Z <- model.matrix(~ Group)
  X <- apply(X_raw, 2, function(v) as.numeric(residuals(lm(v ~ Z - 1))))

  # pure group effect
  G <- matrix(as.integer(Group == "Old"), ncol = 1)
  K_pat <- ncol(X)

  # -----------------------------
  # Shared sample-level latent factors
  # -----------------------------
  if (n_shared_factors < 0) stop("n_shared_factors must be >= 0.")
  L_shared <- if (n_shared_factors == 0) {
    matrix(0, nrow = N, ncol = 0)
  } else {
    matrix(rnorm(N * n_shared_factors, sd = shared_factor_sd),
           nrow = N, ncol = n_shared_factors)
  }

  g <- as.numeric(Group == "Old")
  g <- g - mean(g)              # center
  L_shared[, 1] <- L_shared[, 1] + g

  # -----------------------------
  # Block-wise signal counts (supports scalar / vector / named vector)
  # -----------------------------
  block_names <- names(blocks)
  B <- length(block_names)

  expand_to_blocks <- function(x, what) {
    if (length(x) == 1) {
      out <- rep(as.integer(x), B)
      names(out) <- block_names
      return(out)
    }
    if (!is.null(names(x))) {
      if (!all(block_names %in% names(x))) {
        stop(sprintf("%s: named vector must contain all block names: %s",
                     what, paste(block_names, collapse = ", ")))
      }
      out <- as.integer(x[block_names])
      names(out) <- block_names
      return(out)
    }
    if (length(x) != B) {
      stop(sprintf("%s: must be length 1 or length %d (number of blocks).", what, B))
    }
    out <- as.integer(x)
    names(out) <- block_names
    out
  }

  p_pattern_relevant <- expand_to_blocks(p_pattern_relevant, "p_pattern_relevant")
  p_diff_relevant    <- expand_to_blocks(p_diff_relevant,    "p_diff_relevant")

  for (bn in block_names) {
    p_total <- as.integer(blocks[[bn]])
    if (!is.finite(p_total) || p_total <= 0) stop(sprintf("%s: p_total must be positive.", bn))
    if (p_pattern_relevant[[bn]] < 0) stop(sprintf("%s: p_pattern_relevant must be >= 0.", bn))
    if (p_diff_relevant[[bn]]    < 0) stop(sprintf("%s: p_diff_relevant must be >= 0.", bn))
    if (p_pattern_relevant[[bn]] > p_total) stop(sprintf("%s: p_pattern_relevant > p_total.", bn))
    if (p_diff_relevant[[bn]]    > p_total) stop(sprintf("%s: p_diff_relevant > p_total.", bn))
  }

  # -----------------------------
  # Generate each omics block
  # -----------------------------
  X_blocks <- list()
  truth <- list()
  loadings <- list()
  betas <- list()

  for (bn in block_names) {
    p_total <- as.integer(blocks[[bn]])
    p_pat   <- p_pattern_relevant[[bn]]
    p_diff  <- p_diff_relevant[[bn]]

    # pattern features (block-specific)
    index.pat <- if (p_pat > 0) seq_len(p_pat) else integer(0)
    beta_pat <- matrix(0, nrow = p_total, ncol = K_pat)
    if (p_pat > 0) {
      beta_pat[index.pat, ] <- matrix(
        rnorm(length(index.pat) * K_pat, mean = mean.pat, sd = sd.pat),
        ncol = K_pat
      )
    }

    # differential features (block-specific; allow overlap within block)
    set.seed(999 + match(bn, block_names))
    min.num <- min(p_pat, p_diff)

    if (p_diff == 0) {
      index.diff <- integer(0)
    } else if (p_pat == 0) {
      index.diff <- sample(seq_len(p_total), p_diff)
    } else {
      n_inter <- floor(min.num * percentage_overlap_variables)
      n_inter <- max(min(n_inter, p_pat), 0)

      index.inter <- if (n_inter > 0) sample(index.pat, n_inter) else integer(0)
      remaining <- p_diff - length(index.inter)

      pool <- (1:p_total)[-index.pat]
      if (remaining > length(pool)) stop(sprintf("%s: not enough non-pattern features for diff.only.", bn))

      index.diff.only <- if (remaining > 0) sample(pool, remaining) else integer(0)
      index.diff <- c(index.inter, index.diff.only)
    }

    beta_diff <- rep(0, p_total)
    if (length(index.diff) > 0) {
      beta_diff[index.diff] <- rnorm(length(index.diff), mean = mean.diff, sd = sd.diff)
    }

    # shared latent contribution (block-specific loadings)
    if (n_shared_factors == 0) {
      A_sig <- matrix(0, nrow = 0, ncol = p_total)
      latent_part <- matrix(0, nrow = N, ncol = p_total)
    } else {
      latent_part <- matrix(0, nrow = N, ncol = p_total)

      signal_idx <- sort(unique(c(index.pat, index.diff)))

      if (length(signal_idx) > 0 && n_shared_factors > 0) {
        A_sig <- matrix(
          rnorm(n_shared_factors * length(signal_idx), sd = loading_sd),
          nrow = n_shared_factors
        )
        latent_part[, signal_idx] <- L_shared %*% A_sig
      }

    }

    # signal + noise
    mu <- G %*% t(matrix(beta_diff, ncol = 1)) + X %*% t(beta_pat) + latent_part
    noise <- matrix(rnorm(N * p_total, mean = mean.bg, sd = sd.bg), nrow = N)

    X_blocks[[bn]] <- mu + noise
    colnames(X_blocks[[bn]]) <- paste0("F",1:p_total)
    truth[[bn]] <- list(true.pat = index.pat, true.diff = index.diff)
    loadings[[bn]] <- A_sig
    betas[[bn]] <- list(beta_pat = beta_pat, beta_diff = beta_diff)
  }

  list(
    X_blocks = X_blocks,
    Age = Age,
    Group = Group,
    L_shared = L_shared,
    loadings = loadings,
    truth = truth
  )
}


