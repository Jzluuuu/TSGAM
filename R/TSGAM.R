############################################################################
#' Two-segment generalized additive models for feature screening
#'
#' This function utilizes two-segment generalized additive models to identify
#' variables showing significant differences between young and old segments
#' based on age information.
#'
#'
#'
#' @param X A numeric matrix dimension \eqn{n \times p} representing features
#' measured on samples, or a list including such matrices
#' which multi-omics profiles measured on the same set of samples.
#' @param Age A numeric vector variables of length \eqn{n} containing age
#' information for each sample.
#' @param ncomp Integer, the number of components to include in the DIABLO
#' model.
#' @param keep_rate Numeric, The proportion of features retained by splsda.
#' @param young_cut A numeric value specifying the upper age threshold for the
#' Young segment. Samples with \code{Age <= young_cut} are assigned to the Young
#' segment.
#' @param old_cut A numeric value specifying the lower age threshold for the Old
#' segment. Samples with \code{Age >= old_cut} are assigned to the Old segment.
#' @param p_delta_thr A numeric value specifying the significance threshold
#' for delta.
#' @param adjr2_min A numeric value specifying the minimum adjusted
#' \eqn{R^2} required.
#' @param dev_min A numeric value specifying the minimum proportion of
#' deviance explained required.
#' @param edf_thr A numeric value specifying the threshold for linear or
#' nonlinear.
#' @param feat_k An integer specifying the basis dimension of spline smoothers
#' used in feature-level TS-GAM fitting.
#' @param feat_gamma A numeric value controlling the smoothness penalty
#' inflation in TS-GAM fitting, used to stabilize model estimation.
#' @param feat_nboot An integer specifying the number of bootstrap resamples
#' used to assess the uncertainty and significance of segment-specific effects.
#' @param multi Logical, should be set to \code{TRUE}, if the \code{X} input is
#' a list, the DIABLO is using before the TS-GAM.
#' @param WGCNA Logical, should be set to \code{TRUE}, if the one of \code{X}
#' matrics have high dimension. The WGCNA is using before the feature-level
#' TS-GAM.
#'
#'
#' @return \code{TSGAM} returns a list that contains the following components:
#'
#' \item{res_boot}{A long-format data frame containing bootstrap-refined TS-GAM
#' results for selected candidate features, including bootstrap-based standard
#' errors and p-values for segment-specific effects.}
#' \item{SigPrimary}{A data frame of broadly significant features passing
#' bootstrap-based significance and minimum model fit criteria, representing a
#' conservative set of differential candidates.}
#' \item{SigMarginal}{A data frame of marginal (expression-dominated)
#' differential features, characterized by large segment-specific effects with
#' similar trajectory shapes across age segments.}
#' \item{SigTrajectory}{A data frame of trajectory-differential features
#' exhibiting segment-specific nonlinear age patterns or pronounced differences
#' in trajectory complexity between segments.}

#'
#' @author Jianzhong Lu
#'
#'
#'
#' @examples
#' ## First example
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
#' result <- TSGAM(sim$data,
#'                 sim$Age)
#'
#'
#' @export
TSGAM <- function(X,
                  Age,
                  ncomp = 1,
                  keep_rate = 0.8,
                  young_cut = 36,
                  old_cut = 49,
                  p_delta_thr = 0.05,
                  adjr2_min = 0.30,
                  dev_min   = 0.30,
                  edf_thr = 1.0,
                  feat_k = 6,
                  feat_gamma = 1.2,
                  feat_nboot = 20,
                  multi = FALSE,
                  WGCNA= TRUE){

  # Group info
  ph <- prepare_pheno(Age, young_cut, old_cut)
  Group <- ifelse(ph$Young == 1, "Young", "Old")

  if (class(X)[1] != "list"){
    input_list <- vector(mode='list', 1)
    input_list[[1]] = X
  } else {
    input_list = X
  }

  for (i in 1:length(input_list)){
    idx <- as.integer(sub("^S", "", ph$sample))
    input_list[[i]] <- scale(input_list[[i]][idx,])
  }

  if (multi == TRUE){
    DIABLO_list <- TSGAM_DIABLO(input_list,
                               Age,
                               Group,
                               ncomp,
                               keep_rate)

    for (i in 1:length(input_list)){
      input_list[[i]] <- input_list[[i]][,colnames(input_list[[i]]) %in% rownames(DIABLO_list[[i]])]
    }
  }

  if (WGCNA == TRUE){
      ME_output <- TSGAM_WGCNA(input_list, Age, Group)

      ME_select_feature <- ME_matrix <- vector("list", length(input_list))
      for (i in 1:length(ME_output)){
        if (length(table(ME_output[[i]]$colors)) > 2){
          ME_result <- m_TSGAM(ME_output,ph$Age,Group,
                               k = feat_k, gamma = feat_gamma, bs = "cs")

          ME_select_feature[[i]] <- names(ME_output[[i]][["colors"]])[ME_output[[i]][["colors"]] %in% ME_result[[i]]]
          ME_matrix[[i]] <- input_list[[i]][,colnames(input_list[[i]]) %in% ME_select_feature[[i]]]
        } else {
          ME_matrix[[i]] = input_list[[i]]
        }
      }


  }else {
      ME_matrix = input_list
  }

  light_TSGAM_result <- TSGAM_light(ME_matrix,ph$Age,Group,
                                    k = feat_k, gamma = feat_gamma, bs = "cs")

  light_TSGAM_matrix <- vector(mode='list', length(input_list))
  for (i in 1:length(input_list)){
    light_TSGAM_matrix[[i]] <- input_list[[i]][,colnames(input_list[[i]]) %in% light_TSGAM_result[[i]]]
  }

  boot_result <- TSGAM_boot(light_TSGAM_matrix, ph$Age, Group,
                       k = feat_k, nboot = feat_nboot, edf_thr = edf_thr, gamma = feat_gamma, bs = "cs")

  sig_primary <- sig_diff_expr <- sig_diff_pattern <- vector(mode='list', length(input_list))
  for (i in 1:length(input_list)){
    eff60 <- stats::quantile(abs(boot_result[[i]]$SE_Delta[is.finite(boot_result[[i]]$SE_Delta)]), 0.60, na.rm=TRUE)
    eff80 <- stats::quantile(abs(boot_result[[i]]$SE_Delta[is.finite(boot_result[[i]]$SE_Delta)]), 0.80, na.rm=TRUE)

    sig_primary[[i]] <- boot_result[[i]] |>
      dplyr::filter(
        is.finite(p_delta),
        is.finite(SE_Delta),
        is.finite(adj_r2),
        p_delta < p_delta_thr,
        (adj_r2 >= adjr2_min | dev_expl >= dev_min),
        abs(SE_Delta) >= eff60
      ) |>
      dplyr::arrange(p_delta)


    sig_diff_expr[[i]] <- sig_primary[[i]] |>
      dplyr::filter(abs(SE_Delta) >= eff80) |>
      dplyr::arrange(p_delta)


    sig_diff_pattern[[i]] <- boot_result[[i]] |>
      dplyr::filter(!shape_pattern %in% "other",
             is.finite(adj_r2),
             (adj_r2 >= adjr2_min | dev_expl >= dev_min))
  }

  result <- list(res_boot = boot_result,
                 sig_primary = sig_primary,
                 sig_diff_expr = sig_diff_expr,
                 sig_diff_pattern = sig_diff_pattern)

  return(invisible(result))
}
