# ----------------------------
# phenotype (Age-cut)
# ----------------------------
prepare_pheno <- function(Age, young_cut = 36, old_cut = 49) {

  ph <- data.frame(
    sample   = paste0("S", seq_along(Age)),
    Age      = as.numeric(Age),
    Young = as.integer(as.numeric(Age) <= young_cut),
    Old   = as.integer(as.numeric(Age) >  old_cut),
    stringsAsFactors = FALSE
  )

  ph <- ph[!(ph$Young == 0L & ph$Old == 0L), , drop = FALSE]
  rownames(ph) <- ph$sample

  ph
}



# --------------------------------------
# DIABLO multi-omics integration
# --------------------------------------
#' @import mixOmics
TSGAM_DIABLO <- function(X,
                         Age,
                         Group,
                         ncomp = 1,
                         keep_rate = 0.8){
  ncomp <- ncomp
  AgeClass <- as.factor(Group)

  keepX_list <- vector(mode='list', length(X))
  for (i in 1:length(X)) {
    keepX_list[[i]] <- ceiling(ncol(X[[i]]) * keep_rate)
  }
  names(keepX_list) = names(X)

  design <- matrix(0.1, nrow = length(X), ncol = length(X))
  diag(design) <- 0

  diab <- block.splsda(
    X      = X,
    Y      = AgeClass,
    ncomp  = ncomp,
    keepX  = keepX_list,
    design = design,
    scale  = FALSE
  )

  L_mat <- vector("list", length(X))

  for (i in seq_along(X)) {

    L <- diab$loadings[[i]][, 1:ncomp, drop = FALSE]

    keep <- rowSums(abs(L)) > 0

    L_mat[[i]] <- as.data.frame(L[keep, , drop = FALSE])
    colnames(L_mat[[i]]) <- paste0("comp", seq_len(ncomp))
  }


  return(L_mat)
}


# --------------------------------------
# Network-based module detection (WGCNA)
# --------------------------------------
#' @import WGCNA dynamicTreeCut
TSGAM_WGCNA <- function(L_mat, Age, group){

  res <- vector("list", length(L_mat))


  for (i in 1:length(L_mat)){
    sel_mat <- colnames(L_mat[[i]])[colSums(abs(L_mat[[i]])) > 0]
    datExpr_X <- L_mat[[i]][, sel_mat, drop = FALSE]

    powers <- 1:20
    sft_X  <- suppressWarnings(pickSoftThreshold(datExpr_X, powerVector = powers, verbose = 0))
    beta_X <- sft_X$powerEstimate

    if (is.na(beta_X)) beta_X <- 6

    adj_X     <- adjacency(datExpr_X, power = beta_X, type = "signed")
    TOM_X     <- TOMsimilarity(adj_X, TOMType = "signed")
    dissTOM_X <- 1 - TOM_X

    tree_X <- hclust(as.dist(dissTOM_X), method = "average")

    module_labels_X <- cutreeDynamic(
      dendro            = tree_X,
      distM             = dissTOM_X,
      deepSplit         = 3,
      pamRespectsDendro = FALSE,
      minClusterSize    = 50
    )

    module_colors_X0 <- labels2colors(module_labels_X)
    names(module_colors_X0) <- colnames(datExpr_X)

    if (length(setdiff(unique(module_colors_X0), "grey")) == 0) {
      res[[i]] <- list(MEs = NULL, colors = module_colors_X0)
    } else {
      merge_X <- mergeCloseModules(
        exprData = datExpr_X,
        colors   = module_colors_X0,
        cutHeight = 0.2,
        verbose  = 0
      )

      module_colors_X <- merge_X$colors
      names(module_colors_X) <- colnames(datExpr_X)

      MEs_X <- merge_X$newMEs
      MEs_X <- orderMEs(MEs_X)
      res[[i]] <- list(MEs = MEs_X, colors = module_colors_X)
    }

  }

  return(res)
}


# --------------------------------------
# Moudle-level TS-GAM
# --------------------------------------
#' @import mgcv tibble dplyr stats
m_TSGAM <- function(MEs, Age, Group,
                    k = 6, gamma = 1.0, bs = "cs") {

  df_base <- data.frame(
    Age      = as.numeric(Age),
    YoungSeg = as.integer(Group == "Young"),
    OldSeg   = as.integer(Group == "Old")
  )

  keep_one <- function(y) {
    y <- as.numeric(y)

    df <- df_base
    df$Y <- y

    m <- mgcv::bam(
      Y ~ s(Age, by = YoungSeg, bs = bs, k = k) +
          s(Age, by = OldSeg, bs = bs, k = k),
      data = df, method = "fREML", gamma = gamma, discrete = TRUE
    )


    sm <- summary(m)
    s_tab <- sm$s.table

    edf_y <- if (any(grepl("YoungSeg", rownames(s_tab))))
      s_tab[grep("YoungSeg", rownames(s_tab))[1], "edf"] else NA_real_
    edf_o <- if (any(grepl("OldSeg", rownames(s_tab))))
      s_tab[grep("OldSeg", rownames(s_tab))[1], "edf"] else NA_real_

    pred  <- as.numeric(predict(m, newdata = df))
    Delta <- mean(pred[df$OldSeg == 1], na.rm = TRUE) - mean(pred[df$YoungSeg == 1], na.rm = TRUE)

    m_result <- tibble::tibble(
      Delta = Delta,
      adj_r2 = sm$r.sq,
      dev_expl = sm$dev.expl,
      edf_young = as.numeric(edf_y),
      edf_old   = as.numeric(edf_o)
    )

    m_result
  }


  keep_X <- vector(mode='list', length(MEs))
  for (i in 1:length(MEs)) {
    m_result_all <- matrix(nrow = ncol(MEs[[i]][["MEs"]]),ncol = 5)
    colnames(m_result_all) <- c('Delta','adj_r2','dev_expl','edf_young','edf_old')
    rownames(m_result_all) <- colnames(MEs[[i]][["MEs"]])
    for (j in 1:ncol(MEs[[i]][["MEs"]])){
      m_result_all[j,] <- as.matrix(keep_one(MEs[[i]][["MEs"]][,j]))
    }
    m_result_all <- as.data.frame(m_result_all)
    thr_me <- m_result_all %>%
      dplyr::mutate(edf_diff = abs(edf_young - edf_old)) %>%
      dplyr::summarise(
        expr80     = stats::quantile(abs(Delta), 0.8, na.rm = TRUE),
        edf_diff80 = stats::quantile(edf_diff,   0.8, na.rm = TRUE)
      )

    m_selected <- m_result_all %>%
      dplyr::mutate(edf_diff = abs(edf_young - edf_old)) %>%
      dplyr::filter(
        is.finite(Delta),
        is.finite(edf_diff),
        abs(Delta) >= thr_me$expr80 |
          edf_diff   >= thr_me$edf_diff80
      )
    rownames(m_selected) <- sub("^ME", "", rownames(m_selected))
    keep_X[[i]] <- rownames(m_selected)
  }

  return(keep_X)
}


# ----------------------------
# Feature-level TS-GAM (light)
# ----------------------------
#' @import mgcv tibble dplyr stats
TSGAM_light <- function(Mat, Age, Group,
                        k = 6, gamma = 1.0, bs = "cs", select = FALSE) {

  df_base <- data.frame(
    Age      = as.numeric(Age),
    YoungSeg = as.integer(Group == "Young"),
    OldSeg   = as.integer(Group == "Old")
  )

  keep_light <- function(y) {
    y <- as.numeric(y)

    df <- df_base
    df$Y <- y

    if (sum(is.finite(df$Y)) < 5 || stats::var(df$Y, na.rm = TRUE) == 0) {
      return(tibble(term=NA_character_, edf=NA_real_, pval=NA_real_,
                    Delta=NA_real_, SE_Delta=NA_real_, adj_r2=NA_real_,
                    dev_expl=NA_real_, p_delta=NA_real_, model=NA_character_))
    }

    form <- Y ~ s(Age, by=YoungSeg, bs=bs, k=k) + s(Age, by=OldSeg, bs=bs, k=k)

    m <- mgcv::bam(formula=form, data=df, method="fREML", select=select, gamma=gamma, discrete=TRUE)
    sm <- summary(m)
    st <- sm$s.table

    pred <- as.numeric(predict(m, newdata=df))
    Delta_obs <- mean(pred[df$OldSeg==1L], na.rm=TRUE) - mean(pred[df$YoungSeg==1L], na.rm=TRUE)

    tibble(
      term=rownames(st),
      edf=as.numeric(st[,"edf"]),
      pval=as.numeric(st[,"p-value"]),
      Delta=Delta_obs,
      adj_r2=sm$r.sq,
      dev_expl=sm$dev.expl
    )
  }

  light_result <- keep_X <- vector(mode='list', length(Mat))
  for (i in 1:length(Mat)) {

    light_matrix <- matrix(nrow = 2*ncol(Mat[[i]]),ncol = 6)
    colnames(light_matrix) <- c('term','edf','pval','Delta','adj_r2','dev_expl')

    light_matrix <- as.data.frame(light_matrix)
    for (j in 1:ncol(Mat[[i]])){
      light_matrix[(2*j-1):(2*j),] <-keep_light(Mat[[i]][,j])
    }
    light_matrix$feature <- rep(colnames(Mat[[i]]), each = 2)

    light_result <- light_matrix |>
      dplyr::mutate(
        seg = dplyr::case_when(
          grepl("YoungSeg", term) ~ "young",
          grepl("OldSeg", term)   ~ "old",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(seg)) |>
      dplyr::group_by(feature) |>
      dplyr::summarise(
        Delta      = dplyr::first(Delta),
        adj_r2     = dplyr::first(adj_r2),
        dev_expl   = dplyr::first(dev_expl),
        edf_young  = dplyr::first(edf[seg == "young"]),
        edf_old    = dplyr::first(edf[seg == "old"]),
        pval_young = dplyr::first(pval[seg == "young"]),
        pval_old   = dplyr::first(pval[seg == "old"]),
        .groups = "drop"
      )


    thr_light <- light_result %>%
      dplyr::mutate(edf_diff = abs(edf_young - edf_old)) %>%
      dplyr::summarise(
        expr80     = stats::quantile(abs(Delta), 0.8, na.rm = TRUE),
        edf_diff80 = stats::quantile(edf_diff,   0.8, na.rm = TRUE)
      )

    light_selected <- light_result %>%
      dplyr::mutate(edf_diff = abs(edf_young - edf_old)) %>%
      dplyr::filter(
        is.finite(Delta),
        is.finite(edf_diff),
        abs(Delta) >= thr_light$expr80 |
          edf_diff   >= thr_light$edf_diff80
      )

    keep_X[[i]] <- light_selected$feature
  }

  return(keep_X)

}


# ----------------------------
# Feature-level TS-GAM (boot)
# ----------------------------
#' @import mgcv tibble dplyr stats
TSGAM_boot <- function(Mat, Age, Group,
                             k = 6, nboot = 100, gamma = 1.2, edf_thr = NA,bs = "cs", select = FALSE) {
  df_base <- data.frame(
    Age      = as.numeric(Age),
    YoungSeg = as.integer(Group == "Young"),
    OldSeg   = as.integer(Group == "Old")
  )

  keep_boot <- function(y) {
    y <- as.numeric(y)

    df <- df_base
    df$Y <- y

    if (sum(is.finite(df$Y)) < 5 || stats::var(df$Y, na.rm=TRUE) == 0) {
      return(tibble(term=NA_character_, edf=NA_real_, pval=NA_real_,
                    Delta=NA_real_, SE_Delta=NA_real_, adj_r2=NA_real_,
                    dev_expl=NA_real_, p_delta=NA_real_, model="error"))
    }

    form <- Y ~ s(Age, by=YoungSeg, bs=bs, k=k) + s(Age, by=OldSeg, bs=bs, k=k)
    m <- mgcv::bam(formula=form, data=df, method="fREML", select=select, gamma=gamma, discrete=TRUE)
    sm <- summary(m)

    pred <- as.numeric(predict(m, newdata=df))
    Delta_obs <- mean(pred[df$OldSeg==1L], na.rm=TRUE) - mean(pred[df$YoungSeg==1L], na.rm=TRUE)

    idxY <- which(df$YoungSeg==1L); idxO <- which(df$OldSeg==1L)

    boot_deltas <- replicate(nboot, {
      iY <- sample(idxY, replace=TRUE)
      iO <- sample(idxO, replace=TRUE)
      dfb <- df[c(iY,iO),,drop=FALSE]

      vY <- stats::var(dfb$Y[dfb$YoungSeg==1L], na.rm=TRUE)
      vO <- stats::var(dfb$Y[dfb$OldSeg==1L],   na.rm=TRUE)
      if (!is.finite(vY) || !is.finite(vO) || vY == 0 || vO == 0) return(NA_real_)

      mb <- mgcv::bam(formula=form, data=dfb, method="fREML",
                      select=select, gamma=gamma, discrete=TRUE)
      prb <- as.numeric(predict(mb, newdata=dfb))
      mean(prb[dfb$OldSeg==1L], na.rm=TRUE) - mean(prb[dfb$YoungSeg==1L], na.rm=TRUE)
    })

    boot_deltas <- boot_deltas[is.finite(boot_deltas)]
    SE_Delta <- stats::sd(boot_deltas, na.rm=TRUE)

    p_delta  <- if (is.finite(SE_Delta) && SE_Delta>0) 2*stats::pnorm(-abs(Delta_obs/SE_Delta)) else NA_real_

    st <- sm$s.table
    tibble(
      term=rownames(st), edf=as.numeric(st[,"edf"]), pval=as.numeric(st[,"p-value"]),
      Delta=Delta_obs, SE_Delta=SE_Delta, adj_r2=sm$r.sq, dev_expl=sm$dev.expl,
      p_delta=p_delta
    )
  }


  boot_result <- vector(mode='list', length(Mat))
  for (i in 1:length(Mat)) {

    boot_matrix <- matrix(nrow = 2*ncol(Mat[[i]]),ncol = 8)
    colnames(boot_matrix) <- c('term','edf','pval','Delta','SE_Delta','adj_r2','dev_expl','p_delta')

    boot_matrix <- as.data.frame(boot_matrix)
    for (j in 1:ncol(Mat[[i]])){
      boot_matrix[(2*j-1):(2*j),] <-keep_boot(Mat[[i]][,j])
    }
    boot_matrix$feature <- rep(colnames(Mat[[i]]), each = 2)

    boot_result[[i]] <- boot_matrix |>
      dplyr::mutate(
        seg = dplyr::case_when(
          grepl("YoungSeg", term) ~ "young",
          grepl("OldSeg", term)   ~ "old",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(seg)) |>
      dplyr::group_by(feature) |>
      dplyr::summarise(
        Delta      = Delta[1],
        SE_Delta   = SE_Delta[1],
        p_delta    = p_delta[1],
        adj_r2     = adj_r2[1],
        dev_expl   = dev_expl[1],
        edf_young  = edf[seg=="young"][1],
        edf_old    = edf[seg=="old"][1],
        pval_young = pval[seg=="young"][1],
        pval_old   = pval[seg=="old"][1],
        .groups="drop"
      )

    if(is.na(edf_thr)){edf_thr  <- quantile(c(boot_result[[i]]$edf_young, boot_result[[i]]$edf_old), 0.8, na.rm=TRUE)}
    diff_thr <- quantile(abs(boot_result[[i]]$edf_young - boot_result[[i]]$edf_old), 0.8, na.rm=TRUE)

    boot_result[[i]] <- boot_result[[i]] |>
      mutate(
        young_nonlinear = (edf_young >= edf_thr) & (is.finite(pval_young) & pval_young < 0.05),
        old_nonlinear   = (edf_old   >= edf_thr) & (is.finite(pval_old)   & pval_old   < 0.05),
        edf_diff        = abs(edf_young - edf_old),
        edf_diff_top    = edf_diff >= diff_thr,
        shape_pattern = case_when(
          young_nonlinear & !old_nonlinear ~ "young_nonlinear_old_linear",
          !young_nonlinear & old_nonlinear ~ "young_linear_old_nonlinear",
          TRUE                             ~ "other"
        )
      )
  }

  return(boot_result)
}




