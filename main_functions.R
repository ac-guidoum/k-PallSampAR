###############################################################################
##  PTE_KSAMPLE_parallel() — unique function, equal or unequal n_i 
##
##  n_vec        : vector of sample sizes (n_1,...,n_k) — balanced if rep(n_i,k)
##  rho_1        : ρ₀ under H₀
##  deltat_grid  : grid of δ values to simulate (scalar or vector)
##  M            : number of Monte-Carlo replications
##  ncores       : cores (single cluster for the entire grid)
##  p            : shrinkage weight ∈ (0,1)
##  seed         : reproducible seed
##
##  Returns a tibble (one row per δ): deltat, SRE_RE, …, SRE_JSplus
###############################################################################

PTE_KSAMPLE_parallel <- function(n_vec,
                                  rho_1       = 0.3,
                                  deltat_grid = 0,
                                  M           = 10000,
                                  ncores      = parallel::detectCores() - 1L,
                                  p           = 0.5,
                                  seed        = NULL,
                                  return_z2   = FALSE) {

  k       <- length(n_vec)
  n_total <- sum(n_vec)

  if (k < 2L)           stop("k = length(n_vec) must be >= 2.")
  if (any(n_vec < 3L))  stop("All n_vec elements must be >= 3.")
  if (abs(rho_1) >= 1)  stop("|rho_1| must be < 1.")

  crit_005 <- qchisq(0.95, df = k - 1L)
  crit_01  <- qchisq(0.90, df = k - 1L)

  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  if (!is.null(seed)) doRNG::registerDoRNG(seed)

  rows <- vector("list", length(deltat_grid))

  for (di in seq_along(deltat_grid)) {

    d         <- deltat_grid[di]
    delta_vec <- seq(0, d, length.out = k)
    rho_vec   <- rho_1 + delta_vec / sqrt(n_total)

    if (any(abs(rho_vec) >= 1)) {
      warning(sprintf("deltat = %.4g skipped: |rho_j| >= 1 for some j.", d))
      next
    }

    parallel::clusterExport(
      cl,
      varlist = c("n_vec", "rho_vec", "k", "p", "crit_005", "crit_01"),
      envir   = environment()
    )

    flat <- foreach::foreach(
      m             = 1:M,
      .combine      = c,
      .multicombine = TRUE,
      .packages     = "stats"
    ) %dopar% {

      rho_hat_j <- nu_hat_j <- numeric(k)

      for (j in seq_len(k)) {
        n_j <- n_vec[j]
        x   <- arima.sim(model = list(ar = rho_vec[j]), n = n_j)
        est <- ar.ols(x, aic = FALSE, order.max = 1, demean = TRUE, intercept = FALSE)
        rho_hat_j[j] <- as.numeric(est$ar)
        sigma2        <- est$var.pred
        gamma0        <- mean((x - mean(x))^2)       
        nu_hat_j[j]  <- sigma2 / gamma0
      }

      w         <- n_vec / nu_hat_j
      rhoRE     <- sum(w * rho_hat_j) / sum(w)
      rhoRE_vec <- rep(rhoRE, k)
      dv        <- rho_hat_j - rhoRE_vec
      rhoSE     <- p * rhoRE + (1 - p) * rho_hat_j

      Z2 <- sum(dv^2 * w)
      PE005 <- if (Z2 < crit_005) rhoRE_vec else rho_hat_j
      PE01  <- if (Z2 < crit_01)  rhoRE_vec else rho_hat_j
      SP005 <- if (Z2 < crit_005) rhoSE     else rho_hat_j
      SP01  <- if (Z2 < crit_01)  rhoSE     else rho_hat_j

      if (k >= 4L && Z2 > 0) {
        sh  <- (k - 3L) / Z2
        JS  <- rho_hat_j - sh * dv
        JSp <- rhoRE_vec + max(1 - sh, 0) * dv
      } else {
        JS <- JSp <- rho_hat_j
      }

      c(rho_hat_j, rhoRE_vec, rhoSE, PE005, PE01, SP005, SP01, JS, JSp, Z2)
    }

    nc       <- 9L * k + 1L
    mat      <- matrix(flat, nrow = M, ncol = nc, byrow = TRUE)
    sl       <- function(b) mat[, ((b - 1L) * k + 1L):(b * k), drop = FALSE]
    true_mat <- matrix(rho_vec, nrow = M, ncol = k, byrow = TRUE)
    mse      <- function(b) colMeans((sl(b) - true_mat)^2)
    mue      <- mse(1L)

    out <- tibble::tibble(
      deltat     = d,
      SRE_RE     = mean(mue / mse(2L)),
      SRE_SE     = mean(mue / mse(3L)),
      SRE_PE_005 = mean(mue / mse(4L)),
      SRE_PE_01  = mean(mue / mse(5L)),
      SRE_SP_005 = mean(mue / mse(6L)),
      SRE_SP_01  = mean(mue / mse(7L)),
      SRE_JS     = mean(mue / mse(8L)),
      SRE_JSplus = mean(mue / mse(9L))
    )
    if (return_z2) out$Z2 <- list(mat[, nc])
    rows[[di]] <- out
  }

  dplyr::bind_rows(rows)
}

# Function show resultants

show4 <- function(tab, digits = 4) {
  out <- tab %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("SRE"),
                                ~ formatC(.x, format = "f", digits = digits)))
  print(out, n = Inf)
  invisible(out)
}

# Function plot 

plot_sre <- function(sre_table, n_vec, rho_0 = 0.3) {
  
  k <- length(n_vec)
  ni_val <- if (length(unique(n_vec)) == 1L) {
    unique(n_vec)  
  } else {
    paste0("(", paste(n_vec, collapse = ","), ")")      
  }
  make_label <- function(a) {
    ni_expr <- if (is.numeric(ni_val)) {
      paste0("n[i]==", ni_val)
    } else {
      paste0('n[i]=="', ni_val, '"')
    }
    paste0("rho[0]==", rho_0,
           "~','~", ni_expr,
           "~','~k==", k,
           "~','~alpha==", a)
  }
  
  sre_long <- sre_table %>%
    pivot_longer(cols = -deltat, names_to = "Estimator", values_to = "SRE") %>%
    mutate(
      alpha = case_when(
        Estimator == "SRE_PE_005" ~ 0.05,
        Estimator == "SRE_PE_01"  ~ 0.10,
        Estimator == "SRE_SP_005" ~ 0.05,
        Estimator == "SRE_SP_01"  ~ 0.10,
        TRUE ~ NA_real_
      ),
      Estimator = case_when(
        Estimator == "SRE_RE"     ~ "RE",
        Estimator == "SRE_SE"     ~ "SE",
        Estimator == "SRE_JS"     ~ "JS",
        Estimator == "SRE_JSplus" ~ "JS+",
        Estimator == "SRE_PE_005" ~ "PE",
        Estimator == "SRE_PE_01"  ~ "PE",
        Estimator == "SRE_SP_005" ~ "SP",
        Estimator == "SRE_SP_01"  ~ "SP",
        TRUE ~ Estimator
      )
    )
  
  sre_expanded <- sre_long %>%
    filter(is.na(alpha)) %>%
    dplyr::select(-alpha) %>%
    tidyr::expand_grid(alpha = c(0.05, 0.10)) %>%
    bind_rows(sre_long %>% filter(!is.na(alpha)))
  
  sre_expanded$alpha_label <- factor(
    sre_expanded$alpha,
    levels = c(0.05, 0.10),
    labels = c(make_label(0.05), make_label(0.10))
  )
  
  deltas <- unique(sre_expanded$deltat)
  df_ue  <- data.frame(
    deltat     = rep(deltas, 2),
    Estimator  = "UE",
    SRE        = 1,
    alpha      = rep(c(0.05, 0.10), each = length(deltas)),
    alpha_label = rep(levels(sre_expanded$alpha_label), each = length(deltas))
  )
  
  sre_all <- bind_rows(sre_expanded, df_ue)
  sre_all$Estimator <- factor(
    sre_all$Estimator,
    levels = c("UE", "RE", "SE", "PE", "SP", "JS", "JS+")
  )
  colors_transparent <- c(
    "UE"  = "black",
    "RE"  = alpha("#E41A1C", 0.5),
    "SE"  = alpha("#377EB8", 0.5),
    "PE"  = alpha("#4DAF4A", 0.5),
    "SP"  = alpha("#984EA3", 0.5),
    "JS"  = alpha("#FF7F00", 0.5),
    "JS+" = alpha("#A65628", 0.5)
  )
  
  ggplot(sre_all, aes(x = deltat, y = SRE,
                      color    = Estimator,
                      linetype = Estimator)) +
    geom_line(aes(linewidth = Estimator)) +
    scale_linewidth_manual(
      values = c("UE"=1.1, "RE"=0.9, "SE"=0.9,
                 "PE"=0.9, "SP"=0.9, "JS"=0.9, "JS+"=0.9)
    ) +
    scale_color_manual(values = colors_transparent) +
    scale_linetype_manual(
      values = c("UE"="dashed", "RE"="solid", "SE"="solid",
                 "PE"="solid",  "SP"="solid", "JS"="solid", "JS+"="solid")
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 10),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    scale_y_continuous(breaks = pretty_breaks(n = 10)) +
    facet_wrap(~ alpha_label, labeller = label_parsed) +
    labs(
      x = expression(delta),
      y = "Simulated Relative Efficiency"
    ) +
    guides(
      color     = guide_legend(nrow = 1),
      linetype  = guide_legend(nrow = 1),
      linewidth = guide_legend(nrow = 1)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.75),
      panel.spacing    = unit(0.02, "lines"),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 12, face = "bold"),
      legend.key.width = unit(1.4, "cm"),
      legend.box       = "horizontal",
      legend.box.margin  = margin(t = -10),
      legend.box.spacing = unit(0.2, "cm"),
      axis.title       = element_text(face = "bold"),
      strip.text       = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.75)
    )
}

################################################################################
# Main function to calculate ADR for different estimators
###############################################################################

calculate_ADR <- function(k, Delta, alpha = 0.05, p = 0.5) {
  
  expect_chi2_inv <- function(df, Delta) {
    integrate(function(x) {
      dchisq(x, df = df, ncp = Delta) / x
    }, lower = 1e-10, upper = Inf)$value
  }

  expect_chi4_inv <- function(df, Delta) {
    integrate(function(x) {
      dchisq(x, df = df, ncp = Delta) / (x^2)
    }, lower = 1e-10, upper = Inf)$value
  }

  expect_truncated <- function(df, Delta, c) {
    integrate(function(x) {
      (1 - c/x) * dchisq(x, df = df, ncp = Delta)
    }, lower = 1e-10, upper = c)$value
  }

  expect_truncated_sq <- function(df, Delta, c) {
    integrate(function(x) {
      (1 - c/x)^2 * dchisq(x, df = df, ncp = Delta)
    }, lower = 1e-10, upper = c)$value
  }
  ADR_UE <- k
  ADR_RE <- 1 + Delta
  ADR_SE <- k - p*(2-p)*(k-1) + p^2*Delta
  crit_val <- qchisq(1-alpha, df = k-1)
  G_k1 <- pchisq(crit_val, df = k+1, ncp = Delta)
  G_k3 <- pchisq(crit_val, df = k+3, ncp = Delta)
  ADR_PE <- k - (k-1)*G_k1 + Delta*(2*G_k1 - G_k3)
  ADR_SP <- k - p*(2-p)*(k-1)*G_k1 + p*Delta*(2*G_k1 - (2-p)*G_k3)
  E1 <- expect_chi2_inv(k+1, Delta)
  E2 <- expect_chi4_inv(k+1, Delta)
  E3 <- expect_chi4_inv(k+3, Delta)
  ADR_JS <- k - (k-1)*(k-3)*(2*E1 - (k-3)*E2) + (k-3)*(k+1)*Delta*E3
  c_val <- k-3
  T1 <- expect_truncated_sq(k+1, Delta, c_val)
  T2 <- expect_truncated(k+1, Delta, c_val)
  T3 <- expect_truncated_sq(k+3, Delta, c_val)
  ADR_JSplus <- ADR_JS - (k-1)*T1 + Delta*(2*T2 - T3)
  return(data.frame(
    Delta = Delta,
    UE = ADR_UE,
    RE = ADR_RE,
    SE = ADR_SE,
    PE = ADR_PE,
    SP = ADR_SP,
    JS = ADR_JS,
    JSplus = ADR_JSplus
  ))
}

################################################################################
# Main function to calculate ADQB for different estimators
###############################################################################

calculate_ADQB <- function(k, Delta, alpha = 0.05, p = 0.5) {
  G <- function(x, df, delta) {
    pchisq(x, df, delta)
  }
  E_chi2_inv <- function(df, delta) {
    integrate(function(x) { dchisq(x, df, delta) / x }, lower = 0, upper = Inf)$value
  }
  E_chi2_inv_trunc <- function(df, delta, c) {
    integrate(function(x) { dchisq(x, df, delta) / x }, lower = c, upper = Inf)$value
  }
  q_alpha <- qchisq(1 - alpha, df = k - 1)
  term_G_k1 <- G(q_alpha, k + 1, Delta)
  term_G_k3 <- G(q_alpha, k + 3, Delta)
  term_E_inv <- E_chi2_inv(k + 1, Delta)
  term_E_inv_trunc <- E_chi2_inv_trunc(k + 1, Delta, k - 3)
  term_G_trunc <- G(k - 3, k + 1, Delta)
  return(data.frame(
    Delta = Delta,
    UE = 0,
    RE = Delta,
    SE = p^2 * Delta,
    PE = Delta * term_G_k1^2,
    SP = p^2 * Delta * term_G_k1^2,
    JS = (k - 3)^2 * Delta * term_E_inv^2,
    JSplus = Delta * (term_G_trunc + (k - 3) * term_E_inv_trunc)^2
  ))
}

###############################################################################
# test stationarity
###############################################################################

test_stationarity <- function(x, name = "series", alpha = 0.05){
  x <- stats::na.omit(as.numeric(x))
  
  adf_d <- ur.df(x, type="drift", selectlags="AIC")
  adf_t <- ur.df(x, type="trend", selectlags="AIC")
  
  adf_d_stat <- unname(adf_d@teststat["tau2"])
  adf_t_stat <- unname(adf_t@teststat["tau3"])
  adf_d_cv5  <- unname(adf_d@cval["tau2","5pct"])
  adf_t_cv5  <- unname(adf_t@cval["tau3","5pct"])
  adf_d_reject <- is.finite(adf_d_stat) && adf_d_stat < adf_d_cv5
  adf_t_reject <- is.finite(adf_t_stat) && adf_t_stat < adf_t_cv5
  
  kpss_l <- tseries::kpss.test(x, null="Level")
  kpss_t <- tseries::kpss.test(x, null="Trend")
  kpss_l_reject <- kpss_l$p.value < alpha
  kpss_t_reject <- kpss_t$p.value < alpha
  
  adf_any    <- adf_d_reject || adf_t_reject
  kpss_anyR  <- kpss_l_reject || kpss_t_reject
  kpss_allNR <- (!kpss_l_reject) && (!kpss_t_reject)   
  
  verdict <- if (adf_any && kpss_allNR) {
    "Stationnaire (I(0))"
  } else if (!adf_any && kpss_anyR) {
    "Non stationnaire (I(1))"
  } else {
    "Indécis (proche de la racine unitaire ou rupture possible)"
  }
  
  out <- list(
    name = name,
    ADF_drift = list(stat = adf_d_stat, cv5 = adf_d_cv5, reject = adf_d_reject, lags = adf_d@lags, spec="const"),
    ADF_trend = list(stat = adf_t_stat, cv5 = adf_t_cv5, reject = adf_t_reject, lags = adf_t@lags, spec="const+trend"),
    KPSS_level = list(stat = as.numeric(kpss_l$statistic), p.value = kpss_l$p.value, reject = kpss_l_reject),
    KPSS_trend = list(stat = as.numeric(kpss_t$statistic), p.value = kpss_t$p.value, reject = kpss_t_reject),
    verdict = verdict
  )
  class(out) <- "stationarity_result"
  out
}

print.stationarity_result <- function(x, ...){
  cat(sprintf("=== stationarity de %s ===\n", x$name))
  cat(sprintf("ADF (const):     stat = %.3f, cv5 = %.3f, lags = %d, reject = %s\n",
              x$ADF_drift$stat, x$ADF_drift$cv5, x$ADF_drift$lags, x$ADF_drift$reject))
  cat(sprintf("ADF (trend):     stat = %.3f, cv5 = %.3f, lags = %d, reject = %s\n",
              x$ADF_trend$stat, x$ADF_trend$cv5, x$ADF_trend$lags, x$ADF_trend$reject))
  cat(sprintf("KPSS (level):    stat = %.3f, p = %.4f, reject = %s\n",
              x$KPSS_level$stat, x$KPSS_level$p.value, x$KPSS_level$reject))
  cat(sprintf("KPSS (trend):    stat = %.3f, p = %.4f, reject = %s\n",
              x$KPSS_trend$stat, x$KPSS_trend$p.value, x$KPSS_trend$reject))
  cat(sprintf("\nVerdict: %s\n", x$verdict))
}

# Fonction : teste H0: rho = 1  vs H1: rho < 1 (ADF / DF-GLS) + KPSS
test_root_less_than_one <- function(x, type = c("drift","trend","none"), alpha = 0.05){
  type <- match.arg(type)
  x <- na.omit(as.numeric(x))
  if(length(x) < 10) warning("Très peu d'observations — résultats peu fiables.")
  
  # ADF (ur.df) : null = unit root; si on rejette => preuve rho < 1
  adf <- ur.df(x, type = ifelse(type=="none","none",
                                ifelse(type=="drift","drift","trend")),
               selectlags = "AIC")
  # extraits ADF
  adf_stat <- if (type=="none") adf@teststat["tau1"] else if (type=="drift") adf@teststat["tau2"] else adf@teststat["tau3"]
  adf_cv5  <- if (type=="none") adf@cval["tau1","5pct"] else if (type=="drift") adf@cval["tau2","5pct"] else adf@cval["tau3","5pct"]
  adf_reject <- is.finite(adf_stat) && (adf_stat < adf_cv5)
  
  # DF-GLS (ur.ers) : souvent plus puissant près de rho ~ 1
  # model argument: "constant" or "trend"
  model_ers <- if (type=="trend") "trend" else "constant"
  ers <- tryCatch(ur.ers(x, type = "DF-GLS", model = model_ers, lag.max = floor(12*(length(x)/100)^0.25)),
                  error = function(e) NULL)
  ers_stat <- if (!is.null(ers)) ers@teststat else NA_real_
  ers_reject <- NA
  # ur.ers doesn't always provide cvals in same slot; we will print summary instead of automating cv check
  
  # KPSS (null = stationarity) as a complement
  kpss_level <- suppressWarnings(tryCatch(kpss.test(x, null = "Level"), error=function(e) NULL))
  kpss_reject <- if (!is.null(kpss_level)) (kpss_level$p.value < alpha) else NA
  
  # output
  cat("---- Test résumé ----\n")
  cat(sprintf("ADF (type=%s): stat = %.4f, 5%% cv = %.4f --> reject H0 (unit root)? %s\n",
              type, as.numeric(adf_stat), as.numeric(adf_cv5), ifelse(adf_reject,"YES","NO")))
  if(!is.null(ers)){
    cat("DF-GLS (DF-GLS stat) :\n"); print(summary(ers))
  } else cat("DF-GLS : impossible à calculer (erreur).\n")
  if (!is.null(kpss_level)){
    cat(sprintf("KPSS (level): stat = %.4f, p = %.4f --> reject H0 (stationary)? %s\n",
                as.numeric(kpss_level$statistic), kpss_level$p.value, ifelse(kpss_reject,"YES","NO")))
  } else cat("KPSS : impossible à calculer.\n")
  
  interpretation <- if (adf_reject && !kpss_reject) {
    "Preuve que rho < 1 (série stationnaire) — ADF rejette et KPSS ne rejette pas."
  } else if (!adf_reject && kpss_reject) {
    "Preuve que rho = 1 (non stationnaire) — ADF ne rejette pas et KPSS rejette."
  } else {
    "Indécis : résultats contradictoires ou proches de la frontière (faire DF-GLS, Ng-Perron, test de rupture, ou bootstrap)."
  }
  cat("\nInterpretation: ", interpretation, "\n")
  invisible(list(adf = list(stat = adf_stat, cv5 = adf_cv5, reject = adf_reject),
                 dfgls = ers, kpss = kpss_level, interpretation = interpretation))
}

###############################################################################

fit_chisq <- function(data) {
  nll_chisq_non_central <- function(df, ncp) {
    -sum(dchisq(data, df = df, ncp = ncp, log = TRUE))
  }
  fit <- stats4::mle(nll_chisq_non_central, start = list(df = trunc(var(data)/2), 
                                                         ncp = trunc(mean(data))),
                     lower = list(df = 0.01, ncp = 0), optim = optimx::optimr)
  params <- coef(fit)
  return(params)
}

cal_statics <- function(serie) {
  est      <- ar.ols(serie, aic = FALSE, order.max = 1, demean = FALSE, intercept = FALSE)
  gamma    <- acf(serie, plot = FALSE, lag.max = 0, type = "covariance", demean = TRUE)$acf[1]
  sigma2   <- est$var.pred
  mean_val <- mean(serie)          # renamed: avoids shadowing base::mean
  variance <- var(serie)
  n        <- length(serie)
  nu       <- sigma2 / gamma
  rho      <- est$ar[[1]]
  return(c(n, mean_val, variance, sigma2, gamma, nu, rho))
}

###############################################################################

estimate_estimators <- function(data, shrink_p = 0.5, alpha = c(0.05, 0.10),
                                digits   = 4)
{
  stopifnot(is.data.frame(data))
  ar_stats <- function(x) {
    x <- stats::na.omit(x)
    fit <- stats::ar.ols(x, aic = FALSE, order.max = 1,
                         demean = FALSE, intercept = FALSE)
    sigma2 <- fit$var.pred
    gamma0 <- stats::acf(x, plot = FALSE, type = "covariance",
                         lag.max = 0, demean = TRUE)$acf[1]
    list(rho = as.numeric(fit$ar),
         nu  = sigma2 / gamma0,
         n   = length(x))
  }
  num_df <- data |>
    dplyr::select(where(is.numeric)) |>
    dplyr::select(-year)
  stats  <- purrr::map(num_df, ar_stats)
  rho_hat<- purrr::map_dbl(stats, "rho")
  nu_hat <- purrr::map_dbl(stats, "nu")
  n_vec  <- purrr::map_dbl(stats, "n")
  k      <- length(rho_hat)
  w      <- n_vec / nu_hat
  rho_RE <- sum(w * rho_hat) / sum(w)
  rho_SE <- shrink_p * rho_RE + (1 - shrink_p) * rho_hat
  Z2     <- sum(w * (rho_hat - rho_RE)^2)
  crit   <- stats::qchisq(1 - alpha, df = k - 1)
  names(crit) <- paste0("crit_", sprintf("%03d", alpha*100))
  rho_PE <- lapply(crit, \(cval) if (Z2 < cval) rep(rho_RE, k) else rho_hat)
  rho_SP <- lapply(crit, \(cval) if (Z2 < cval) rho_SE       else rho_hat)
  shrink      <- (k - 3) / Z2
  rho_JS      <- rho_hat - shrink * (rho_hat - rho_RE)
  shrink_pos  <- max(1 - (k - 3) / Z2, 0)
  rho_JSpos   <- rho_RE + shrink_pos * (rho_hat - rho_RE)
  est_tab <- dplyr::tibble(
    Country = names(rho_hat),
    UE      = rho_hat,
    RE      = rep(rho_RE, k),
    SE      = rho_SE,
    PE_005  = rho_PE[[1]],
    PE_010  = rho_PE[[2]],
    SP_005  = rho_SP[[1]],
    SP_010  = rho_SP[[2]],
    JS      = rho_JS,
    JSpos   = rho_JSpos
  ) 
  pretty_tab <- est_tab |>
    dplyr::mutate(
      dplyr::across(where(is.numeric),
                    \(x) sprintf(paste0("%.", digits, "f"), x))
    )
  
  list(table = pretty_tab,       
       Z2    = round(Z2, digits))
}
