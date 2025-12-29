#' Confidence Intervals for AUC Under MAR Verification
#'
#' Computes point estimates and confidence intervals for the AUC of a continuous
#' test when disease verification is missing at random (MAR). The function returns
#' four estimates simultaneously, obtained using the bias-corrected estimators FI, MSI,
#' IPW, and SPE proposed by Alonzo and Pepe (2005).
#'
#' @param T Test results; a positive numeric vector.
#' @param D Verified disease status; a logical vector with possible missing values.
#' @param A Covariate; a positive numeric vector. Only one covariate is allowed.
#' @param alpha Significance level for the confidence interval. Default is 0.05.
#' @param search_step Step size used in root searching. Default is 0.01.
#' @param tol Tolerance used in root searching. Default is 1e-5.
#' @param precision Precision parameter used in the regression model. Default is 1e-4.
#' @param n.boot Number of bootstrap replicates. Default is 1000.
#' @param plot Logical; if `TRUE` (default) a density plot is produced.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{n.total}}{Total number of subjects.}
#'   \item{\code{n.case}}{Number of verified diseased subjects.}
#'   \item{\code{n.control}}{Number of verified non-diseased subjects.}
#'   \item{\code{p.missing}}{Proportion of missing verification.}
#'   \item{\code{pt.est}}{Point estimates of AUC.}
#'   \item{\code{BC.intervals}}{Bootstrap classic (BC) confidence intervals.}
#'   \item{\code{BP.intervals}}{Bootstrap percentile (BP) confidence intervals.}
#'   \item{\code{HEL1.intervals}}{Hybrid empirical likelihood confidence intervals, type I.}
#'   \item{\code{HEL2.intervals}}{Hybrid empirical likelihood confidence intervals, type II.}
#' }
#'
#' @references
#' Alonzo, T. A. and Pepe, M. S. (2005). Assessing accuracy of a continuous
#' screening test in the presence of verification bias. *Journal of the Royal
#' Statistical Society: Series C (Applied Statistics)*.
#'
#' Wang, S., Shi, S., and Qin, G. (2025). Empirical likelihood inference for the area
#' under the ROC curve with verification-biased data. Manuscript under peer review.
#'
#' @details
#' Bootstrap and hybrid empirical likelihood confidence intervals for AUC under
#' verification bias are computed.
#'
#' The disease model \eqn{\rho} is estimated using a probit regression model
#' linear in \eqn{T} and \eqn{A} based on verified subjects, given by
#'
#' \deqn{
#'   \rho_i = P(D_i = 1 \mid T_i, A_i)
#'   = \Phi(\alpha + \beta T_i + \gamma A_i),
#'   \quad i = 1, \ldots, n.
#' }
#'
#' where \eqn{\Phi} denotes the standard normal cumulative distribution function.
#'
#' The verification model is estimated using a logit regression model
#' linear in \eqn{T} and \eqn{A} based on all subjects, given by
#'
#' \deqn{
#'   \operatorname{logit}(\pi_i)
#'   = \log\!\left( \frac{\pi_i}{1 - \pi_i} \right)
#'   = \alpha + \beta T_i + \gamma A_i,
#'   \quad i = 1, \ldots, n,
#' }
#'
#' where \eqn{\pi_i = P(V_i = 1 \mid T_i, A_i)}.
#'
#' The function may also produce a density plot of the test measurements when `plot = TRUE`.
#'
#' @examples
#' set.seed(123)
#' T <- abs(rnorm(100))
#' A <- abs(rnorm(100))
#' D <- as.logical(T + A > stats::quantile(T + A, 0.8))
#' D[sample(100, 30)] <- NA
#' auc.ci.mar(T, D, A, n.boot = 20, plot = FALSE)

#' @export
auc.ci.mar <- function(T,
                       D,
                       A,
                       alpha = 0.05,
                       search_step = 0.01,
                       tol = 1e-5,
                       precision = 1e-4,
                       n.boot = 1000,
                       plot = TRUE) {
  bc.scale = function(x) {
    if (any(is.na(x)))
      stop("Input contains NA values.")
    if (any(x <= 0))
      stop("Input not positive.")
    bc <- MASS::boxcox(x ~ 1, plotit = FALSE)
    # The optimal lambda is the value on the x-axis corresponding to the peak
    # of the log-likelihood curve on the y-axis.
    lambda_optimal <- bc$x[which.max(bc$y)]

    # --- Apply Transformation ---
    # A lambda of 0 corresponds to a natural log transformation.
    if (lambda_optimal == 0) {
      transformed_x <- log(x)
    } else {
      transformed_x <- (x^lambda_optimal - 1) / lambda_optimal
    }
    return(drop(scale(transformed_x)))
  }
  #Function used to calculate weights given sample
  get_weights = function(T, D, A, V, n) {
    # Fit Disease model
    disease = stats::glm(
      D[V == 1] ~ T[V == 1] + A[V == 1],
      family = stats::binomial(link = 'probit'),
      control = list(maxit = 100)
    )
    beta = disease$coefficients
    rho = drop(stats::pnorm(cbind(rep(1, n), T, A) %*% beta))
    rho[rho > 1 - precision] = 1 - precision
    rho[rho < precision] = precision
    # Fit Verification model
    Verif = stats::glm(
      V ~ T + A,
      family = stats::binomial(link = 'logit'),
      control = list(maxit = 100)
    )
    betaV = Verif$coefficients
    pi_hat = drop(1 / (1 + exp(-(
      cbind(1, T, A) %*% betaV
    ))))
    pi_hat[pi_hat > 1 - precision] = 1 - precision
    pi_hat[pi_hat < precision] = precision
    # Calculate weights
    wd.fi  = rho
    wd.msi = V * D + (1 - V) * rho
    wd.ipw = V * D / pi_hat
    wd.spe = (V * D - (V - pi_hat) * rho) / pi_hat
    wn.fi  = 1 - rho
    wn.msi = V * (1 - D) + (1 - V) * (1 - rho)
    wn.ipw = V * (1 - D) / pi_hat
    wn.spe = (V * (1 - D) - (V - pi_hat) * (1 - rho)) / pi_hat
    return(list(
      wd  = rbind(wd.fi, wd.msi, wd.ipw, wd.spe),
      wn  = rbind(wn.fi, wn.msi, wn.ipw, wn.spe),
      rho = rho,
      pi  = pi_hat
    ))
  }
  #Calculate bootstrap Copies of auc
  BSP = function(dummy) {
    boots = sample(1:n, n, replace = TRUE)
    Tb = T[boots]
    Db = D[boots]
    Ab = A[boots]
    Vb = V[boots]
    VD = sum(Vb * Db)
    VH = sum(Vb * (1 - Db))
    if (VD >= 1 && VH >= 1)
    {
      wtsb = get_weights(Tb, Db, Ab, Vb, n)
      #calculate delta
      TT = matrix(rep(Tb, n), n)
      Ib = (TT >= t(TT)) + 0
      ghatb = (wtsb$wd %*% Ib) / apply(wtsb$wd, 1, sum)
      deltab = (ghatb %*% t(wtsb$wn)) / apply(wtsb$wn, 1, sum)
      #deltab gives all 16 combs, we only need 4 diagonal elements
      HEL_Ub = drop(wtsb$wn %*% Ib / apply(wtsb$wn, 1, sum))
      HEL_Lb = wtsb$wd * (1 - HEL_Ub - deltahat) / apply(wtsb$wd, 1, sum)
      HEL_rb = apply(HEL_Lb, 1, function(row) {
        tryCatch(
          emplik::el.test(row, mu = 0)$'-2LLR',
          error = function(e)
            NA
        )
      })
      return(c(diag(deltab), unname(HEL_rb)))
    }
  }

  HEL_LLR = function(delta) {
    #calculate delta and U for HEL
    #HEL_U = drop(wts$wn[i, ] %*% It / sum(wts$wn[i, ]))
    HEL_L = wts$wd[i, ] * (1 - HEL_U[i, ] - delta) / sum(wts$wd[i, ])
    if (all(HEL_L == 0)) {
      HEL_L = 0.01 + wts$wd[i, ] * (1 - HEL_U[i, ] - delta) / sum(wts$wd[i, ])
    }
    HEL_r = emplik::el.test(HEL_L, mu = 0)$'-2LLR'
    return(HEL_r)
  }
  HEL1_LLR = function(delta) {
    HEL_LLR(delta) - BSQ[i]
  }
  HEL2_LLR = function(delta) {
    HEL_LLR(delta) * (7 / 9)^3 / BSmedian[i] - stats::qchisq(1 - alpha, 1)
  }
  Find_zero = function(f) {
    est_zero = which(sapply(seq(0, 1, search_step), f) < 0)
    if (length(est_zero) == 0) {
      LB = NA
      UB = NA
    } else {
      if (f(0) < 0) {
        LB = 0
      } else {
        LB = stats::uniroot(
          f,
          lower = search_step * (est_zero[1] - 2),
          upper = search_step * (est_zero[1] - 1),
          tol = tol
        )$root
      }
      if (f(1) < 0) {
        UB = 1
      } else {
        UB = stats::uniroot(
          f,
          lower = search_step * (max(est_zero) - 1),
          upper = search_step * max(est_zero),
          tol = tol
        )$root
      }
    }
    return(c(LB, UB))
  }
  V = !is.na(D)
  if (sum(D, na.rm = TRUE) <= 1 ||
      sum(!D, na.rm = TRUE) <= 1)
    stop("Not enough verified subjects (<=1)")
  if (sum(D, na.rm = TRUE) <= 10 ||
      sum(!D, na.rm = TRUE) <= 10)
    warning("Not enough verified subjects (<=10)")
  verified = data.frame(Tv = T[V == 1], Dv = factor(D[V == 1]))
  Density.plot <- ggplot2::ggplot(verified, ggplot2::aes(x = Tv, linetype = Dv)) +
    ggplot2::geom_density(color = "black", linewidth = 1) +
    ggplot2::scale_linetype_manual(
      values = c('TRUE' = "solid", 'FALSE' = "dashed"),
      name = "Group",
      labels = c("TRUE" = "Diseased", "FALSE" = "Not Diseased")
    ) +
    ggplot2::labs(x = 'Measurement', y = "Density") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.width = grid::unit(1.5, "cm"),
      # Wider lines
      legend.key.height = grid::unit(0, "cm"),
      # Slim height
      #legend.spacing.x = grid::unit(0.1, "cm"),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      )
    )
  if (plot)
    plot(Density.plot)
  n = length(T)
  za = stats::qnorm(1 - alpha / 2)
  T = bc.scale(T)
  A = bc.scale(A)
  D = as.numeric(D)
  D[is.na(D)] = 0
  V = as.numeric(V)
  d = (-1)^(pROC::roc(response = verified$Dv, predictor = verified$Tv)$d == ">")
  T = d * T
  #Calculate weights
  wts = get_weights(T, D, A, V, n)
  sum.wd = apply(wts$wd, 1, sum)
  sum.wn = apply(wts$wn, 1, sum)
  TT = matrix(rep(T, n), n)
  It = (TT >= t(TT)) + 0
  ghat = (wts$wd %*% It) / apply(wts$wd, 1, sum)
  deltahat = diag((ghat %*% t(wts$wn)) / apply(wts$wn, 1, sum))
  HEL_U = (wts$wn %*% It) / apply(wts$wn, 1, sum)
  #Bootstrap Intervals
  BSP.est = do.call(cbind, Filter(Negate(is.null), lapply(1:n.boot, FUN = BSP)))
  BSP.delta = BSP.est[1:4, ]
  BSP.r = BSP.est[5:8, ]
  BSP.delta[BSP.delta < 0] = 0
  BSP.delta[BSP.delta > 1] = 1
  BPL = apply(BSP.delta,
              1,
              stats::quantile,
              probs = alpha / 2,
              na.rm = TRUE)
  BPU = apply(BSP.delta,
              1,
              stats::quantile,
              probs = 1 - alpha / 2,
              na.rm = TRUE)
  BCmean = apply(BSP.delta, 1, mean, na.rm = TRUE)
  BCvar = apply(BSP.delta, 1, stats::var, na.rm = TRUE)
  BCL = BCmean - za * sqrt(BCvar)
  BCU = BCmean + za * sqrt(BCvar)
  BCL = pmin(pmax(BCL, 0), 1)
  BCU = pmin(pmax(BCU, 0), 1)
  BSQ = apply(BSP.r,
              1,
              stats::quantile,
              probs = 1 - alpha,
              na.rm = TRUE)
  BSmedian = apply(BSP.r, 1, stats::median, na.rm = TRUE)
  HEL1_CI = NULL
  for (i in 1:4) {
    HEL1_CI = rbind(HEL1_CI, Find_zero(HEL1_LLR))
  }
  HEL2_CI = NULL
  for (i in 1:4) {
    HEL2_CI = rbind(HEL2_CI, Find_zero(HEL2_LLR))
  }
  BC = cbind(BCL, BCU)
  BP = cbind(BPL, BPU)
  HELI = HEL1_CI
  HELII = HEL2_CI
  names(deltahat) = c("FI", "MSI", "IPW", "SPE")
  rownames(BP) = c("FI", "MSI", "IPW", "SPE")
  rownames(BC) = c("FI", "MSI", "IPW", "SPE")
  rownames(HELI) = c("FI", "MSI", "IPW", "SPE")
  rownames(HELII) = c("FI", "MSI", "IPW", "SPE")
  colnames(BP) = c("Lower", "Upper")
  colnames(BC) = c("Lower", "Upper")
  colnames(HELI) = c("Lower", "Upper")
  colnames(HELII) = c("Lower", "Upper")
  return(
    list(
      n.total = n,
      n.case = sum(V * D, na.rm = TRUE),
      n.control = sum(V * (1 - D), na.rm = TRUE),
      p.missing = 1 - sum(V) / n,
      pt.est = deltahat,
      BC.intervals = BC,
      BP.intervals = BP,
      HEL1.intervals = HELI,
      HEL2.intervals = HELII
    )
  )
}

#' Confidence Intervals for Sensitivity at Fixed Level of Specificity Under MAR Verification
#'
#' Computes point estimates and confidence intervals for sensitivity of a continuous test
#' at a fixed level of specificity when disease verification is missing at random (MAR).
#' The function returns four estimates simultaneously, obtained using the bias-corrected estimators
#' FI, MSI, IPW, and SPE proposed by Alonzo and Pepe (2005).
#'
#' @param T Test results; a positive numeric vector.
#' @param D Verified disease status; a logical vector with possible missing values.
#' @param A Covariate; a positive numeric vector. Only one covariate is allowed.
#' @param p Target specificity level; a number between 0 and 1.
#' @param alpha Significance level for the confidence interval. Default is 0.05.
#' @param search_step Step size used in root searching. Default is 0.01.
#' @param tol Tolerance used in root searching. Default is 1e-5.
#' @param precision Precision parameter used in the regression model. Default is 1e-4.
#' @param n.boot Number of bootstrap replicates. Default is 1000.
#' @param plot Logical; if `TRUE` (default) a density plot is produced.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{n.total}}{Total number of subjects.}
#'   \item{\code{n.case}}{Number of verified diseased subjects.}
#'   \item{\code{n.control}}{Number of verified non-diseased subjects.}
#'   \item{\code{p.missing}}{Proportion of missing verification.}
#'   \item{\code{pt.est}}{Point estimates of sensitivity at specificity `p`.}
#'   \item{\code{pt.est.ac}}{Point estimates of sensitivity at specificity `p` using the Agresti--Coull method.}
#'   \item{\code{AC.intervals}}{Agresti--Coull-based confidence intervals.}
#'   \item{\code{WS.intervals}}{Wilson score-based confidence intervals.}
#'   \item{\code{BTI.intervals}}{Bootstrap confidence intervals, type I.}
#'   \item{\code{BTII.intervals}}{Bootstrap confidence intervals, type II.}
#'   \item{\code{HEL1.intervals}}{Hybrid empirical likelihood confidence intervals, type I.}
#'   \item{\code{HEL2.intervals}}{Hybrid empirical likelihood confidence intervals, type II.}
#' }
#'
#' @references
#' Alonzo, T. A. and Pepe, M. S. (2005). Assessing accuracy of a continuous
#' screening test in the presence of verification bias. *Journal of the Royal
#' Statistical Society: Series C (Applied Statistics)*.
#'
#' Wang, S., Shi, S., and Qin, G. (2025). Empirical likelihood-based confidence
#' intervals for sensitivity of a continuous test at a fixed level of specificity
#' with verification bias. Manuscript under peer review.
#'
#' @details
#' The function targets sensitivity evaluated at specificity level `p` (i.e.,
#' sensitivity at the threshold achieving specificity `p`). Bootstrap and
#' hybrid empirical likelihood confidence intervals are computed as returned in the list.
#'
#' The disease model \eqn{\rho} is estimated using a probit regression model
#' linear in \eqn{T} and \eqn{A} based on verified subjects, given by
#'
#' \deqn{
#'   \rho_i = P(D_i = 1 \mid T_i, A_i)
#'   = \Phi(\alpha + \beta T_i + \gamma A_i),
#'   \quad i = 1, \ldots, n.
#' }
#'
#' where \eqn{\Phi} denotes the standard normal cumulative distribution function.
#'
#' The verification model is estimated using a logit regression model
#' linear in \eqn{T} and \eqn{A} based on all subjects, given by
#'
#' \deqn{
#'   \operatorname{logit}(\pi_i)
#'   = \log\!\left( \frac{\pi_i}{1 - \pi_i} \right)
#'   = \alpha + \beta T_i + \gamma A_i,
#'   \quad i = 1, \ldots, n,
#' }
#'
#' where \eqn{\pi_i = P(V_i = 1 \mid T_i, A_i)}.
#'
#' The function may also produce a density plot of the test measurements when `plot = TRUE`.
#'
#' @examples
#' set.seed(123)
#' T <- abs(rnorm(100))
#' A <- abs(rnorm(100))
#' D <- as.logical(T + A > stats::quantile(T + A, 0.8))
#' D[sample(100, 30)] <- NA
#' sen.ci.mar(T, D, A, p = 0.8, n.boot = 20, plot = FALSE)
#'
#' @export
sen.ci.mar <- function(T,
                       D,
                       A,
                       p,
                       alpha = 0.05,
                       search_step = 0.01,
                       tol = 1e-5,
                       precision = 1e-4,
                       n.boot = 1000,
                       plot = TRUE) {
  bc.scale = function(x) {
    if (any(is.na(x)))
      stop("Input contains NA values.")
    if (any(x <= 0))
      stop("Input not positive.")
    bc <- MASS::boxcox(x ~ 1, plotit = FALSE)
    # The optimal lambda is the value on the x-axis corresponding to the peak
    # of the log-likelihood curve on the y-axis.
    lambda_optimal <- bc$x[which.max(bc$y)]

    # --- Apply Transformation ---
    # A lambda of 0 corresponds to a natural log transformation.
    if (lambda_optimal == 0) {
      transformed_x <- log(x)
    } else {
      transformed_x <- (x^lambda_optimal - 1) / lambda_optimal
    }
    return(drop(scale(transformed_x)))
  }
  #Function used to calculate weights given sample
  get_weights = function(T, D, A, V, n) {
    # Fit Disease model
    disease = stats::glm(
      D[V == 1] ~ T[V == 1] + A[V == 1],
      family = stats::binomial(link = 'probit'),
      control = list(maxit = 100)
    )
    beta = disease$coefficients
    rho = drop(stats::pnorm(cbind(rep(1, n), T, A) %*% beta))
    rho[rho > 1 - precision] = 1 - precision
    rho[rho < precision] = precision
    # Fit Verification model
    Verif = stats::glm(
      V ~ T + A,
      family = stats::binomial(link = 'logit'),
      control = list(maxit = 100)
    )
    betaV = Verif$coefficients
    pi_hat = drop(1 / (1 + exp(-(
      cbind(1, T, A) %*% betaV
    ))))
    pi_hat[pi_hat > 1 - precision] = 1 - precision
    pi_hat[pi_hat < precision] = precision
    # Calculate weights
    wd.fi  = rho
    wd.msi = V * D + (1 - V) * rho
    wd.ipw = V * D / pi_hat
    wd.spe = (V * D - (V - pi_hat) * rho) / pi_hat
    wn.fi  = 1 - rho
    wn.msi = V * (1 - D) + (1 - V) * (1 - rho)
    wn.ipw = V * (1 - D) / pi_hat
    wn.spe = (V * (1 - D) - (V - pi_hat) * (1 - rho)) / pi_hat
    return(list(
      wd  = rbind(wd.fi, wd.msi, wd.ipw, wd.spe),
      wn  = rbind(wn.fi, wn.msi, wn.ipw, wn.spe),
      rho = rho,
      pi  = pi_hat
    ))
  }
  #Function used to Calculate Bootstrap Copies of Jtilde
  BSP = function(dummy) {
    boots = sample(1:n, n, replace = TRUE)
    Tb = T[boots]
    Db = D[boots]
    Ab = A[boots]
    Vb = V[boots]
    if (sum(Vb * Db) >= 1 && sum(Vb * (1 - Db)) >= 1)
    {
      wtsb = get_weights(Tb, Db, Ab, Vb, n)
      sum.bwd = apply(wtsb$wd, 1, sum)
      sum.bwn = apply(wtsb$wn, 1, sum)
      #Find J and cutoff
      Tbsort = sort(Tb)
      cbtest = Tbsort[1:n - 1] + (Tbsort[2:n] - Tbsort[1:n - 1]) / 2
      jb = NULL
      for (c in cbtest) {
        tb_temp_2 = (Tb < c) + 0
        #Spe
        speb = drop(tb_temp_2 %*% t(wtsb$wn)) / sum.bwn
        speb = pmin(pmax(speb, 0), 1)
        jb = rbind(jb, speb)
      }
      cpb = cbtest[apply(jb <= p, 2, which.min)]
      #cp = ctest[apply(abs(jm - p), 2, which.min)]
      #Sen
      compb = outer(
        Tb,
        cpb,
        FUN = function(x, y)
          x >= y
      )
      #senb = diag((wtsb$wd %*% compb) / apply(wtsb$wd, 1, sum))
      senb = diag((wtsb$wd %*% compb + 0.5 * za^2) / (apply(wtsb$wd, 1, sum) + za^2))
      TTb = matrix(rep(Tb, n), n)
      Ib = (TTb >= t(TTb)) + 0
      HEL_Ub = drop(wtsb$wn %*% Ib / apply(wtsb$wn, 1, sum))
      HEL_Lb = wtsb$wd * ((HEL_Ub <= 1 - p) - senhat) / apply(wtsb$wd, 1, sum)
      HEL_rb = apply(HEL_Lb, 1, function(row) {
        tryCatch(
          emplik::el.test(row, mu = 0)$'-2LLR',
          error = function(e)
            NA
        )
      })
      return(c(senb, unname(HEL_rb)))
    }
  }
  HEL_LLR = function(delta) {
    #calculate delta and U for HEL
    #HEL_U = drop(wts$wn[i, ] %*% It / sum(wts$wn[i, ]))
    HEL_L = wts$wd[i, ] * ((HEL_U[i, ] <= 1 - p) - as.vector(delta)) / sum(wts$wd[i, ])
    if (all(HEL_L == 0)) {
      HEL_L = wts$wd[i, ] * ((HEL_U[i, ] <= 1 - p) - as.vector(delta) + 0.01) / sum(wts$wd[i, ])
    }
    HEL_r = emplik::el.test(HEL_L, mu = 0)$'-2LLR'
    return(HEL_r)
  }
  HEL1_LLR = function(delta) {
    HEL_LLR(delta) - BSQ[i]
  }
  HEL2_LLR = function(delta) {
    HEL_LLR(delta) * (7 / 9)^3 / BSmedian[i] - stats::qchisq(1 - alpha, 1)
  }
  Find_zero = function(f) {
    est_zero = which(sapply(seq(0, 1, search_step), f) < 0)
    if (length(est_zero) == 0) {
      LB = NA
      UB = NA
    } else {
      if (f(0) < 0) {
        LB = 0
      } else {
        LB = stats::uniroot(
          f,
          lower = search_step * (est_zero[1] - 2),
          upper = search_step * (est_zero[1] - 1),
          tol = tol
        )$root
      }
      if (f(1) < 0) {
        UB = 1
      } else {
        UB = stats::uniroot(
          f,
          lower = search_step * (max(est_zero) - 1),
          upper = search_step * max(est_zero),
          tol = tol
        )$root
      }
    }
    return(c(LB, UB))
  }
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0 || p >= 1) {
    stop("`p` must be a single numeric value strictly between 0 and 1.", call. = FALSE)
  }
  V = !is.na(D)
  if (sum(D, na.rm = TRUE) <= 1 ||
      sum(!D, na.rm = TRUE) <= 1)
    stop("Not enough verified subjects (<=1)")
  if (sum(D, na.rm = TRUE) <= 10 ||
      sum(!D, na.rm = TRUE) <= 10)
    warning("Not enough verified subjects (<=10)")
  verified = data.frame(Tv = T[V == 1], Dv = factor(D[V == 1]))
  Density.plot = ggplot2::ggplot(verified, ggplot2::aes(x = Tv, linetype = Dv)) +
    ggplot2::geom_density(color = "black", linewidth = 1) +
    ggplot2::scale_linetype_manual(
      values = c('TRUE' = "solid", 'FALSE' = "dashed"),
      name = "Group",
      labels = c("TRUE" = "Diseased", "FALSE" = "Not Diseased")
    ) +
    ggplot2::labs(x = 'Measurement', y = "Density") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.width = grid::unit(1.5, "cm"),
      # Wider lines
      legend.key.height = grid::unit(0, "cm"),
      # Slim height
      #legend.spacing.x = grid::unit(0.1, "cm"),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      )
    )
  if (plot)
    plot(Density.plot)
  n = length(T)
  za = stats::qnorm(1 - alpha / 2)
  T = bc.scale(T)
  A = bc.scale(A)
  D = as.numeric(D)
  D[is.na(D)] = 0
  V = as.numeric(V)
  d = (-1)^(pROC::roc(response = verified$Dv, predictor = verified$Tv)$d == ">")
  T = d * T
  #Calculate weights
  wts = get_weights(T, D, A, V, n)
  sum.wd = apply(wts$wd, 1, sum)
  sum.wn = apply(wts$wn, 1, sum)
  TT = matrix(rep(T, n), n)
  It = (TT >= t(TT)) + 0
  HEL_U = (wts$wn %*% It) / apply(wts$wn, 1, sum)
  #Find J and cutoff
  Tsort = sort(T)
  ctest = Tsort[1:n - 1] + (Tsort[2:n] - Tsort[1:n - 1]) / 2
  jm = NULL
  for (c in ctest) {
    t_temp = (T < c) + 0
    #Spe
    spe = drop(t_temp %*% t(wts$wn)) / sum.wn
    spe = pmin(pmax(spe, 0), 1)
    jm = rbind(jm, spe)
  }
  cp = ctest[apply(jm <= p, 2, which.min)]
  comp = outer(
    T,
    cp,
    FUN = function(x, y)
      x >= y
  )
  #Find senhat and sentilde
  senhat = diag((wts$wd %*% comp) / sum.wd)
  sent = diag((wts$wd %*% comp + 0.5 * za^2) / (sum.wd + za^2))
  sent = pmin(pmax(sent, 0), 1)
  #AC&WS Interval
  AC.var = sent * (1 - sent) / (sum.wd + za^2)
  AC.lb = sent - za * sqrt(AC.var)
  AC.ub = sent + za * sqrt(AC.var)
  AC.lb = pmin(pmax(AC.lb, 0), 1)
  AC.ub = pmin(pmax(AC.ub, 0), 1)
  AC.l = AC.ub - AC.lb
  WS.lb = (sent + 0.5 * za^2 / sum.wd - za * sqrt(sent * (1 - sent) / sum.wd +
                                                    0.25 * za^2 / sum.wd^2)) / (1 + za^2 / sum.wd)
  WS.ub = (sent + 0.5 * za^2 / sum.wd + za * sqrt(sent * (1 - sent) / sum.wd +
                                                    0.25 * za^2 / sum.wd^2)) / (1 + za^2 / sum.wd)
  WS.lb = pmin(pmax(WS.lb, 0), 1)
  WS.ub = pmin(pmax(WS.ub, 0), 1)
  WS.l = WS.ub - WS.lb
  #Bootstrap Intervals
  BSP.est = do.call(cbind, Filter(Negate(is.null), lapply(1:n.boot, FUN = BSP)))
  BSP.sen = BSP.est[1:4, ]
  BSP.r = BSP.est[5:8, ]
  #BPL = apply(BSP.est, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
  #BPU = apply(BSP.est, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  BCmean = apply(BSP.sen, 1, mean, na.rm = TRUE)
  BCvar = apply(BSP.sen, 1, stats::var, na.rm = TRUE)
  BCL1 = sent - za * sqrt(BCvar)
  BCU1 = sent + za * sqrt(BCvar)
  BCL1 = pmin(pmax(BCL1, 0), 1)
  BCU1 = pmin(pmax(BCU1, 0), 1)
  BCL2 = BCmean - za * sqrt(BCvar)
  BCU2 = BCmean + za * sqrt(BCvar)
  BCL2 = pmin(pmax(BCL2, 0), 1)
  BCU2 = pmin(pmax(BCU2, 0), 1)
  LBC1 = BCU1 - BCL1
  LBC2 = BCU2 - BCL2
  #LBP = BPU - BPL
  BSQ = apply(BSP.r,
              1,
              stats::quantile,
              probs = 1 - alpha,
              na.rm = TRUE)
  BSmedian = apply(BSP.r, 1, stats::median, na.rm = TRUE)
  HEL1_CI = NULL
  for (i in 1:4) {
    HEL1_CI = rbind(HEL1_CI, Find_zero(HEL1_LLR))
  }
  HEL2_CI = NULL
  for (i in 1:4) {
    HEL2_CI = rbind(HEL2_CI, Find_zero(HEL2_LLR))
  }
  AC = cbind(AC.lb, AC.ub)
  WS = cbind(WS.lb, WS.ub)
  BTI = cbind(BCL1, BCU1)
  BTII = cbind(BCL2, BCU2)
  HELI = HEL1_CI
  HELII = HEL2_CI
  names(senhat) = c("FI", "MSI", "IPW", "SPE")
  names(sent) = c("FI", "MSI", "IPW", "SPE")
  rownames(AC) = c("FI", "MSI", "IPW", "SPE")
  rownames(WS) = c("FI", "MSI", "IPW", "SPE")
  rownames(BTI) = c("FI", "MSI", "IPW", "SPE")
  rownames(BTII) = c("FI", "MSI", "IPW", "SPE")
  rownames(HELI) = c("FI", "MSI", "IPW", "SPE")
  rownames(HELII) = c("FI", "MSI", "IPW", "SPE")
  colnames(AC) = c("Lower", "Upper")
  colnames(WS) = c("Lower", "Upper")
  colnames(BTI) = c("Lower", "Upper")
  colnames(BTII) = c("Lower", "Upper")
  colnames(HELI) = c("Lower", "Upper")
  colnames(HELII) = c("Lower", "Upper")
  return(
    list(
      n.total = n,
      n.case = sum(V * D, na.rm = TRUE),
      n.control = sum(V * (1 - D), na.rm = TRUE),
      p.missing = 1 - sum(V) / n,
      pt.est = senhat,
      pt.est.ac = sent,
      AC.intervals = AC,
      WS.intervals = WS,
      BTI.intervals = BTI,
      BTII.intervals = BTII,
      HEL1.intervals = HELI,
      HEL2.intervals = HELII
    )
  )
}

#' Confidence Intervals for Youden Index Under MAR Verification
#'
#' Computes point estimates and confidence intervals for maximum Youden index of a
#' continuous test when disease verification is missing at random (MAR). The function
#' returns four estimates simultaneously, obtained using the bias-corrected estimators
#' FI, MSI, IPW, and SPE proposed by Alonzo and Pepe (2005).
#'
#' @param T Test results; a positive numeric vector.
#' @param D Verified disease status; a logical vector with possible missing values.
#' @param A Covariate; a positive numeric vector. Only one covariate is allowed.
#' @param alpha Significance level for the confidence interval. Default is 0.05.
#' @param precision Precision parameter used in the regression model. Default is 1e-4.
#' @param n.boot Number of bootstrap replicates. Default is 1000.
#' @param plot Logical; if `TRUE` (default) a density plot is produced.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{n.total}}{Total number of subjects.}
#'   \item{\code{n.case}}{Number of verified diseased subjects.}
#'   \item{\code{n.control}}{Number of verified non-diseased subjects.}
#'   \item{\code{p.missing}}{Proportion of missing verification.}
#'   \item{\code{pt.est}}{Point estimates of the maximum Youden index.}
#'   \item{\code{pt.est.ac}}{Point estimates of the maximum Youden index using the Agresti--Coull method.}
#'   \item{\code{optimal.cutoff}}{Optimal cutoff point of test results that maximizes the Youden index.}
#'   \item{\code{Wald.intervals}}{Wald confidence intervals.}
#'   \item{\code{BCI.intervals}}{Bootstrap classic confidence intervals, type I.}
#'   \item{\code{BCII.intervals}}{Bootstrap classic confidence intervals, type II.}
#'   \item{\code{BPac.intervals}}{Bootstrap percentile confidence intervals.}
#'   \item{\code{MOVERac.intervals}}{MOVER confidence intervals using the Agresti--Coull method.}
#'   \item{\code{MOVERws.intervals}}{MOVER confidence intervals using the Wilson score method.}
#' }
#'
#' @references
#' Alonzo, T. A. and Pepe, M. S. (2005). Assessing accuracy of a continuous
#' screening test in the presence of verification bias. *Journal of the Royal
#' Statistical Society: Series C (Applied Statistics)*.
#'
#' Wang, S., Shi, S., and Qin, G. (2025). Interval estimation for the Youden index
#' of a continuous diagnostic test with verification biased data. *Statistical Methods in Medical Research*.
#'
#' @details
#' Bootstrap and MOVER-based confidence intervals are computed for the maximum Youden index.
#'
#' The disease model \eqn{\rho} is estimated using a probit regression model
#' linear in \eqn{T} and \eqn{A} based on verified subjects, given by
#'
#' \deqn{
#'   \rho_i = P(D_i = 1 \mid T_i, A_i)
#'   = \Phi(\alpha + \beta T_i + \gamma A_i),
#'   \quad i = 1, \ldots, n.
#' }
#'
#' where \eqn{\Phi} denotes the standard normal cumulative distribution function.
#'
#' The verification model is estimated using a logit regression model
#' linear in \eqn{T} and \eqn{A} based on all subjects, given by
#'
#' \deqn{
#'   \operatorname{logit}(\pi_i)
#'   = \log\!\left( \frac{\pi_i}{1 - \pi_i} \right)
#'   = \alpha + \beta T_i + \gamma A_i,
#'   \quad i = 1, \ldots, n,
#' }
#'
#' where \eqn{\pi_i = P(V_i = 1 \mid T_i, A_i)}.
#'
#' The function may also produce a density plot of the test measurements when `plot = TRUE`.
#'
#' @examples
#' set.seed(123)
#' T <- abs(rnorm(100))
#' A <- abs(rnorm(100))
#' D <- as.logical(T + A > stats::quantile(T + A, 0.8))
#' D[sample(100, 30)] <- NA
#' yi.ci.mar(T, D, A, n.boot = 20, plot = FALSE)
#'
#' @export
yi.ci.mar <- function(T,
                      D,
                      A,
                      alpha = 0.05,
                      precision = 1e-4,
                      n.boot = 1000,
                      plot = TRUE) {
  bc.scale = function(x) {
    if (any(is.na(x)))
      stop("Input contains NA values.")
    if (any(x <= 0))
      stop("Input not positive.")
    bc <- MASS::boxcox(x ~ 1, plotit = FALSE)
    # The optimal lambda is the value on the x-axis corresponding to the peak
    # of the log-likelihood curve on the y-axis.
    lambda_optimal <- bc$x[which.max(bc$y)]

    # --- Apply Transformation ---
    # A lambda of 0 corresponds to a natural log transformation.
    if (lambda_optimal == 0) {
      transformed_x <- log(x)
    } else {
      transformed_x <- (x^lambda_optimal - 1) / lambda_optimal
    }
    return(drop(scale(transformed_x)))
  }
  #Function used to calculate weights given sample
  get_weights = function(T, D, A, V, n) {
    # Fit Disease model
    disease = stats::glm(
      D[V == 1] ~ T[V == 1] + A[V == 1],
      family = stats::binomial(link = 'probit'),
      control = list(maxit = 100)
    )
    beta = disease$coefficients
    rho = drop(stats::pnorm(cbind(rep(1, n), T, A) %*% beta))
    rho[rho > 1 - precision] = 1 - precision
    rho[rho < precision] = precision
    # Fit Verification model
    Verif = stats::glm(
      V ~ T + A,
      family = stats::binomial(link = 'logit'),
      control = list(maxit = 100)
    )
    betaV = Verif$coefficients
    pi_hat = drop(1 / (1 + exp(-(
      cbind(1, T, A) %*% betaV
    ))))
    pi_hat[pi_hat > 1 - precision] = 1 - precision
    pi_hat[pi_hat < precision] = precision
    # Calculate weights
    wd.fi  = rho
    wd.msi = V * D + (1 - V) * rho
    wd.ipw = V * D / pi_hat
    wd.spe = (V * D - (V - pi_hat) * rho) / pi_hat
    wn.fi  = 1 - rho
    wn.msi = V * (1 - D) + (1 - V) * (1 - rho)
    wn.ipw = V * (1 - D) / pi_hat
    wn.spe = (V * (1 - D) - (V - pi_hat) * (1 - rho)) / pi_hat
    return(list(
      wd  = rbind(wd.fi, wd.msi, wd.ipw, wd.spe),
      wn  = rbind(wn.fi, wn.msi, wn.ipw, wn.spe),
      rho = rho,
      pi  = pi_hat
    ))
  }
  #Function used to Calculate Bootstrap Copies of Jtilde
  BSP = function(dummy) {
    boots = sample(1:n, n, replace = TRUE)
    Tb = T[boots]
    Db = D[boots]
    Ab = A[boots]
    Vb = V[boots]
    if (sum(Vb * Db) >= 1 && sum(Vb * (1 - Db)) >= 1)
    {
      wtsb = get_weights(Tb, Db, Ab, Vb, n)
      sum.bwd = apply(wtsb$wd, 1, sum)
      sum.bwn = apply(wtsb$wn, 1, sum)
      Gtildeb = (diag(wtsb$wd %*% drop(outer(Tb, cp, '<'))) + 0.5 * za^2) / (sum.bwd + za^2)
      Ftildeb = (diag(wtsb$wn %*% drop(outer(Tb, cp, '<'))) + 0.5 * za^2) / (sum.bwn + za^2)
      Gtildeb = pmin(pmax(Gtildeb, 0), 1)
      Ftildeb = pmin(pmax(Ftildeb, 0), 1)
      jtildeb = Ftildeb - Gtildeb
      jtildeb = pmin(pmax(jtildeb, 0), 1)
      return(jtildeb)
    }
  }
  V = !is.na(D)
  if (sum(D, na.rm = TRUE) <= 1 ||
      sum(!D, na.rm = TRUE) <= 1)
    stop("Not enough verified subjects (<=1)")
  if (sum(D, na.rm = TRUE) <= 10 ||
      sum(!D, na.rm = TRUE) <= 10)
    warning("Not enough verified subjects (<=10)")
  verified = data.frame(Tv = T[V == 1], Dv = factor(D[V == 1]))
  Density.plot = ggplot2::ggplot(verified, ggplot2::aes(x = Tv, linetype = Dv)) +
    ggplot2::geom_density(color = "black", linewidth = 1) +
    ggplot2::scale_linetype_manual(
      values = c('TRUE' = "solid", 'FALSE' = "dashed"),
      name = "Group",
      labels = c("TRUE" = "Diseased", "FALSE" = "Not Diseased")
    ) +
    ggplot2::labs(x = 'Measurement', y = "Density") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.width = grid::unit(1.5, "cm"),
      # Wider lines
      legend.key.height = grid::unit(0, "cm"),
      # Slim height
      #legend.spacing.x = grid::unit(0.1, "cm"),
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      )
    )
  if (plot)
    plot(Density.plot)
  n = length(T)
  za = stats::qnorm(1 - alpha / 2)
  T.raw = T
  T = bc.scale(T)
  A = bc.scale(A)
  D = as.numeric(D)
  D[is.na(D)] = 0
  V = as.numeric(V)
  d = (-1)^(pROC::roc(response = verified$Dv, predictor = verified$Tv)$d == ">")
  T = d * T
  #Calculate weights
  wts = get_weights(T, D, A, V, n)
  sum.wd = apply(wts$wd, 1, sum)
  sum.wn = apply(wts$wn, 1, sum)
  #Search for j & cutoff point
  Tsort = sort(T)
  ctest = Tsort[1:n - 1] + (Tsort[2:n] - Tsort[1:n - 1]) / 2
  sen.m = NULL
  spe.m = NULL
  for (k in 1:length(ctest)) {
    c = ctest[k]
    t_temp_1 = (T >= c) + 0
    t_temp_2 = (T < c) + 0
    # Spe
    spe = drop(t_temp_2 %*% t(wts$wn)) / sum.wn
    spe = pmin(pmax(spe, 0), 1)
    # Sen
    sen = drop(t_temp_1 %*% t(wts$wd)) / sum.wd
    sen = pmin(pmax(sen, 0), 1)
    sen.m = rbind(sen.m, sen)
    spe.m = rbind(spe.m, spe)
  }
  jm = apply(sen.m + spe.m - 1, 2, max)
  loc = apply(sen.m + spe.m - 1, 2, which.max)
  cp = ctest[loc]
  cp.raw = apply(rbind(sort(T.raw, decreasing = (d == -1))[loc], sort(T.raw, decreasing = (d == -1))[loc + 1]), 2, mean)
  #Find Fhat & Ghat
  Fhat = spe.m[cbind(loc, 1:4)]
  Ghat = 1 - sen.m[cbind(loc, 1:4)]
  Wald.var = Ghat * (1 - Ghat) / sum.wd + Fhat * (1 - Fhat) / sum.wn
  Wald.lb = jm - za * sqrt(Wald.var)
  Wald.ub = jm + za * sqrt(Wald.var)
  Wald.lb = pmin(pmax(Wald.lb, 0), 1)
  Wald.ub = pmin(pmax(Wald.ub, 0), 1)
  #MOVER methods
  Gtilde = (diag(wts$wd %*% drop(outer(T, cp, '<'))) + 0.5 * za^2) / (sum.wd + za^2)
  Ftilde = (diag(wts$wn %*% drop(outer(T, cp, '<'))) + 0.5 * za^2) / (sum.wn + za^2)
  Gtilde = pmin(pmax(Gtilde, 0), 1)
  Ftilde = pmin(pmax(Ftilde, 0), 1)
  LGac = Gtilde - za * sqrt(Gtilde * (1 - Gtilde) / (sum.wd + za^2))
  UGac = Gtilde + za * sqrt(Gtilde * (1 - Gtilde) / (sum.wd + za^2))
  LFac = Ftilde - za * sqrt(Ftilde * (1 - Ftilde) / (sum.wn + za^2))
  UFac = Ftilde + za * sqrt(Ftilde * (1 - Ftilde) / (sum.wn + za^2))
  LGws = (Ghat + 0.5 * za^2 / sum.wd - za * sqrt(Ghat * (1 - Ghat) / sum.wd +
                                                   0.25 * za^2 / sum.wd^2)) / (1 + za^2 / sum.wd)
  UGws = (Ghat + 0.5 * za^2 / sum.wd + za * sqrt(Ghat * (1 - Ghat) / sum.wd +
                                                   0.25 * za^2 / sum.wd^2)) / (1 + za^2 / sum.wd)
  LFws = (Fhat + 0.5 * za^2 / sum.wn - za * sqrt(Fhat * (1 - Fhat) / sum.wn +
                                                   0.25 * za^2 / sum.wn^2)) / (1 + za^2 / sum.wn)
  UFws = (Fhat + 0.5 * za^2 / sum.wn + za * sqrt(Fhat * (1 - Fhat) / sum.wn +
                                                   0.25 * za^2 / sum.wn^2)) / (1 + za^2 / sum.wn)
  Mover.ac.lb = jm - sqrt((Fhat - LFac)^2 + (UGac - Ghat)^2)
  Mover.ac.ub = jm + sqrt((UFac - Fhat)^2 + (Ghat - LGac)^2)
  Mover.ws.lb = jm - sqrt((Fhat - LFws)^2 + (UGws - Ghat)^2)
  Mover.ws.ub = jm + sqrt((UFws - Fhat)^2 + (Ghat - LGws)^2)
  Mover.ac.lb = pmin(pmax(Mover.ac.lb, 0), 1)
  Mover.ac.ub = pmin(pmax(Mover.ac.ub, 0), 1)
  Mover.ws.lb = pmin(pmax(Mover.ws.lb, 0), 1)
  Mover.ws.ub = pmin(pmax(Mover.ws.ub, 0), 1)
  jhat = Fhat - Ghat
  jtilde = Ftilde - Gtilde
  jtilde = pmin(pmax(jtilde, 0), 1)
  #Bootstrap Intervals
  BSP.est = do.call(cbind, Filter(Negate(is.null), lapply(1:n.boot, FUN = BSP)))
  BPL = apply(BSP.est,
              1,
              stats::quantile,
              probs = alpha / 2,
              na.rm = TRUE)
  BPU = apply(BSP.est,
              1,
              stats::quantile,
              probs = 1 - alpha / 2,
              na.rm = TRUE)
  BCmean = apply(BSP.est, 1, mean, na.rm = TRUE)
  BCvar = apply(BSP.est, 1, stats::var, na.rm = TRUE)
  BCL1 = jtilde - za * sqrt(BCvar)
  BCU1 = jtilde + za * sqrt(BCvar)
  BCL1 = pmin(pmax(BCL1, 0), 1)
  BCU1 = pmin(pmax(BCU1, 0), 1)
  BCL2 = BCmean - za * sqrt(BCvar)
  BCU2 = BCmean + za * sqrt(BCvar)
  BCL2 = pmin(pmax(BCL2, 0), 1)
  BCU2 = pmin(pmax(BCU2, 0), 1)
  Wald = cbind(Wald.lb, Wald.ub)
  BP = cbind(BPL, BPU)
  BCI = cbind(BCL1, BCU1)
  BCII = cbind(BCL2, BCU2)
  MOVER.ac = cbind(Mover.ac.lb, Mover.ac.ub)
  MOVER.ws = cbind(Mover.ws.lb, Mover.ws.ub)
  names(cp.raw) = c("FI", "MSI", "IPW", "SPE")
  names(jhat) = c("FI", "MSI", "IPW", "SPE")
  names(jtilde) = c("FI", "MSI", "IPW", "SPE")
  rownames(Wald) = c("FI", "MSI", "IPW", "SPE")
  rownames(BP) = c("FI", "MSI", "IPW", "SPE")
  rownames(BCI) = c("FI", "MSI", "IPW", "SPE")
  rownames(BCII) = c("FI", "MSI", "IPW", "SPE")
  rownames(MOVER.ac) = c("FI", "MSI", "IPW", "SPE")
  rownames(MOVER.ws) = c("FI", "MSI", "IPW", "SPE")
  colnames(Wald) = c("Lower", "Upper")
  colnames(BP) = c("Lower", "Upper")
  colnames(BCI) = c("Lower", "Upper")
  colnames(BCII) = c("Lower", "Upper")
  colnames(MOVER.ac) = c("Lower", "Upper")
  colnames(MOVER.ws) = c("Lower", "Upper")
  return(
    list(
      n.total = n,
      n.case = sum(V * D, na.rm = TRUE),
      n.control = sum(V * (1 - D), na.rm = TRUE),
      p.missing = 1 - sum(V) / n,
      pt.est = jhat,
      pt.est.ac = jtilde,
      optimal.cutoff = cp.raw,
      Wald.intervals = Wald,
      BPac.intervals = BP,
      BCI.intervals = BCI,
      BCII.intervals = BCII,
      MOVERac.intervals = MOVER.ac,
      MOVERws.intervals = MOVER.ws
    )
  )
}
