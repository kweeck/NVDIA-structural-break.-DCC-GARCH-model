# ==============================================================================
# Semiconductor Sector Dynamics: Pre/Post-2015 NVIDIA Structural Break
# Methods: VAR, Granger Causality, IRF, DCC-GJR-GARCH, Rolling Correlation
# Data: Weekly log-returns of NVDA, AMD, INTC, SOXX (2001-2026)
# ==============================================================================

library("aTSA")
library("TSA")
library("forecast")
library("urca")
library("tseries")
library("lmtest")
library("rugarch")
library("fGarch")
library("strucchange")
library("vars")
library("Matrix")
library("matlib")
library("portes")
library("quantmod")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rmgarch")
library("FinTS")
library("zoo")

# ==============================================================================
# DATA: Download and compute log-returns
# ==============================================================================

getSymbols("NVDA", from = "2001-07-09", to = Sys.Date(), periodicity = "weekly")
getSymbols("INTC", from = "2001-07-09", to = Sys.Date(), periodicity = "weekly")
getSymbols("AMD",  from = "2001-07-09", to = Sys.Date(), periodicity = "weekly")
getSymbols("SOXX", from = "2001-07-09", to = Sys.Date(), periodicity = "weekly")

nvd_stat <- diff(log(NVDA$NVDA.Close))
amd_stat <- diff(log(AMD$AMD.Close))
int_stat  <- diff(log(INTC$INTC.Close))
sox_stat  <- diff(log(SOXX$SOXX.Close))

# Split at structural break date
break_date <- "2015-07-27"

nvd_pre  <- na.omit(nvd_stat[paste0("/", break_date)])
nvd_post <- na.omit(nvd_stat[paste0(break_date, "/")])
amd_pre  <- na.omit(amd_stat[paste0("/", break_date)])
amd_post <- na.omit(amd_stat[paste0(break_date, "/")])
int_pre  <- na.omit(int_stat[paste0("/",  break_date)])
int_post <- na.omit(int_stat[paste0(break_date, "/")])
sox_pre  <- na.omit(sox_stat[paste0("/",  break_date)])
sox_post <- na.omit(sox_stat[paste0(break_date, "/")])

# Data frames for VAR / DCC
df_pre <- data.frame(
  NVDA = as.numeric(nvd_pre),
  AMD  = as.numeric(amd_pre),
  INTC = as.numeric(int_pre),
  SOXX = as.numeric(sox_pre)
)

df_post <- data.frame(
  NVDA = as.numeric(nvd_post),
  AMD  = as.numeric(amd_post),
  INTC = as.numeric(int_post),
  SOXX = as.numeric(sox_post)
)

# ==============================================================================
# BLOCK 1: VAR MODEL
# ==============================================================================

# Lag selection — all four criteria (AIC, HQ, SC, FPE) unanimously select VAR(1)
VARselect(df_pre,  lag.max = 20, type = "const")
VARselect(df_post, lag.max = 20, type = "const")

# Estimate VAR(1) for each sub-period
var_pre  <- VAR(df_pre,  p = 1, type = "const")
var_post <- VAR(df_post, p = 1, type = "const")

coef(var_pre)
coef(var_post)

# Residual diagnostics: Portmanteau tests for no cross-correlation
Hosking(var_pre,  lags = 2)
LiMcLeod(var_pre, lags = 2)

Hosking(var_post,  lags = 2)
LiMcLeod(var_post, lags = 2)

# CUSUM stability test: confirms coefficients are stable within each sub-period
plot(stability(var_pre))
plot(stability(var_post))

# ==============================================================================
# BLOCK 2: GRANGER CAUSALITY
# ==============================================================================

# --- System-level causality (each variable vs all others) ---

# Pre-2015
causality(var_pre, cause = "NVDA")$Granger
causality(var_pre, cause = "AMD")$Granger
causality(var_pre, cause = "INTC")$Granger
causality(var_pre, cause = "SOXX")$Granger

# Post-2015
causality(var_post, cause = "NVDA")$Granger
causality(var_post, cause = "AMD")$Granger
causality(var_post, cause = "INTC")$Granger
causality(var_post, cause = "SOXX")$Granger

# --- Pairwise Granger tests: NVDA -> each asset ---
grangertest(as.numeric(amd_pre)  ~ as.numeric(nvd_pre),  order = 1)
grangertest(as.numeric(int_pre)  ~ as.numeric(nvd_pre),  order = 1)
grangertest(as.numeric(sox_pre)  ~ as.numeric(nvd_pre),  order = 1)

grangertest(as.numeric(amd_post) ~ as.numeric(nvd_post), order = 1)
grangertest(as.numeric(int_post) ~ as.numeric(nvd_post), order = 1)
grangertest(as.numeric(sox_post) ~ as.numeric(nvd_post), order = 1)

# --- Reverse direction: does the sector predict NVDA? ---
grangertest(as.numeric(nvd_pre)  ~ as.numeric(sox_pre),  order = 1)
grangertest(as.numeric(nvd_post) ~ as.numeric(sox_post), order = 1)

# --- AMD -> NVDA ---
grangertest(as.numeric(nvd_pre)  ~ as.numeric(amd_pre),  order = 1)
grangertest(as.numeric(nvd_post) ~ as.numeric(amd_post), order = 1)

# ==============================================================================
# BLOCK 3: IMPULSE RESPONSE FUNCTIONS (reduced-form VAR)
# ==============================================================================

extract_irf <- function(irf_obj, period) {
  n <- nrow(irf_obj$irf[[1]])
  data.frame(
    horizon = 0:(n - 1),
    irf     = as.numeric(irf_obj$irf[[1]]),
    lower   = as.numeric(irf_obj$Lower[[1]]),
    upper   = as.numeric(irf_obj$Upper[[1]]),
    period  = period
  )
}

# NVDA shock on other variables — pre-break
irf_pre_nvda_amd  <- irf(var_pre, impulse = "NVDA", response = "AMD",  n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)
irf_pre_nvda_intc <- irf(var_pre, impulse = "NVDA", response = "INTC", n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)
irf_pre_nvda_soxx <- irf(var_pre, impulse = "NVDA", response = "SOXX", n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)

# NVDA shock on other variables — post-break
irf_post_nvda_amd  <- irf(var_post, impulse = "NVDA", response = "AMD",  n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)
irf_post_nvda_intc <- irf(var_post, impulse = "NVDA", response = "INTC", n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)
irf_post_nvda_soxx <- irf(var_post, impulse = "NVDA", response = "SOXX", n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)

# SOXX shock on NVDA — pre and post
irf_pre_soxx_nvda  <- irf(var_pre,  impulse = "SOXX", response = "NVDA", n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)
irf_post_soxx_nvda <- irf(var_post, impulse = "SOXX", response = "NVDA", n.ahead = 12, ortho = FALSE, boot = TRUE, ci = 0.95)

# --- Plot: NVDA -> SOXX ---
df_nvda_soxx <- rbind(
  extract_irf(irf_pre_nvda_soxx,  "Before 2015"),
  extract_irf(irf_post_nvda_soxx, "After 2015")
)

ggplot(df_nvda_soxx, aes(x = horizon, y = irf, color = period, fill = period)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Before 2015" = "#2196F3", "After 2015" = "#E53935")) +
  scale_fill_manual(values  = c("Before 2015" = "#2196F3", "After 2015" = "#E53935")) +
  labs(
    title  = "IRF: NVDA shock \u2192 SOXX response",
    x      = "Weeks ahead",
    y      = "Response",
    color  = "Period",
    fill   = "Period"
  ) +
  theme_minimal()

# --- Plot: SOXX -> NVDA ---
df_soxx_nvda <- rbind(
  extract_irf(irf_pre_soxx_nvda,  "Before 2015"),
  extract_irf(irf_post_soxx_nvda, "After 2015")
)

ggplot(df_soxx_nvda, aes(x = horizon, y = irf, color = period, fill = period)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Before 2015" = "#2196F3", "After 2015" = "#E53935")) +
  scale_fill_manual(values  = c("Before 2015" = "#2196F3", "After 2015" = "#E53935")) +
  labs(
    title = "IRF: SOXX shock \u2192 NVDA response",
    x     = "Weeks ahead",
    y     = "Response",
    color = "Period",
    fill  = "Period"
  ) +
  theme_minimal()

irf_pre_nvda_soxx
irf_post_nvda_soxx

# ==============================================================================
# BLOCK 4: SIGN-RESTRICTED SVAR + FEVD
# ==============================================================================

# --- 4a. Rubio-Ramirez algorithm (random rotations) ---
# Identification via sign restrictions avoids the arbitrary ordering
# imposed by Cholesky decomposition.

sign_restrictions_svar <- function(var_model, sign_mat, n_draws = 10000, n_ahead = 5) {

  Sigma <- summary(var_model)$covres
  P     <- t(chol(Sigma))
  n     <- nrow(Sigma)
  A     <- Acoef(var_model)[[1]]

  accepted_irf <- list()
  n_accept     <- 0

  for (i in seq_len(n_draws)) {

    # Random orthogonal matrix via QR decomposition
    Z <- matrix(rnorm(n * n), n, n)
    Q <- qr.Q(qr(Z))

    # Align column signs so the diagonal of B is positive
    B_cand <- P %*% Q
    signs  <- sign(diag(B_cand))
    signs[signs == 0] <- 1
    Q      <- Q %*% diag(signs)
    B_cand <- P %*% Q

    rownames(B_cand) <- colnames(var_model$y)

    # Check sign restrictions
    check <- TRUE
    for (r in seq_len(nrow(sign_mat))) {
      for (c in seq_len(ncol(sign_mat))) {
        if (!is.na(sign_mat[r, c])) {
          if (sign(B_cand[r, c]) != sign_mat[r, c]) {
            check <- FALSE
            break
          }
        }
      }
      if (!check) break
    }

    if (check) {
      n_accept <- n_accept + 1

      # Compute full IRF for accepted rotation
      irf_mat      <- array(0, dim = c(n, n, n_ahead + 1))
      irf_mat[,,1] <- B_cand
      Phi          <- diag(n)

      for (h in seq_len(n_ahead)) {
        Phi            <- A %*% Phi
        irf_mat[,,h+1] <- Phi %*% B_cand
      }

      accepted_irf[[n_accept]] <- irf_mat
    }
  }

  cat("Accepted:", n_accept, "of", n_draws,
      "(", round(n_accept / n_draws * 100, 2), "%)\n")

  if (n_accept == 0) stop("No rotations accepted — relax sign restrictions")

  irf_array <- simplify2array(accepted_irf)

  irf_mean  <- apply(irf_array, 1:3, mean)
  irf_lower <- apply(irf_array, 1:3, quantile, probs = 0.16)
  irf_upper <- apply(irf_array, 1:3, quantile, probs = 0.84)

  vnames <- colnames(var_model$y)
  dimnames(irf_mean)  <- list(vnames, vnames, 0:n_ahead)
  dimnames(irf_lower) <- list(vnames, vnames, 0:n_ahead)
  dimnames(irf_upper) <- list(vnames, vnames, 0:n_ahead)

  return(list(
    irf_mean  = irf_mean,
    irf_lower = irf_lower,
    irf_upper = irf_upper,
    n_accept  = n_accept,
    n_draws   = n_draws
  ))
}

# --- 4b. Sign restrictions matrix ---
# Variable order: NVDA, AMD, INTC, SOXX
# Rows = responses, columns = structural shocks
# +1 = positive impact required, NA = unrestricted
#
# Economic rationale:
#   NVDA_shock:   positive for NVDA and AMD (AI convergence), unrestricted for INTC/SOXX
#   AMD_shock:    positive for AMD, unrestricted otherwise
#   INTC_shock:   positive for INTC, unrestricted otherwise
#   Sector_shock: positive for all (broad semiconductor shock)

sign_mat <- matrix(c(
#  NVDA  AMD   INTC  Sector
    +1,   NA,   NA,   +1,   # NVDA
    +1,   +1,   NA,   +1,   # AMD
    NA,   NA,   +1,   +1,   # INTC
    NA,   NA,   NA,   +1    # SOXX
), nrow = 4, ncol = 4, byrow = TRUE)

rownames(sign_mat) <- c("NVDA", "AMD", "INTC", "SOXX")
colnames(sign_mat) <- c("NVDA_shock", "AMD_shock", "INTC_shock", "Sector_shock")

# --- 4c. Estimate sign-restricted SVAR ---
set.seed(42)
svar_sign_pre  <- sign_restrictions_svar(var_pre,  sign_mat, n_draws = 10000)
svar_sign_post <- sign_restrictions_svar(var_post, sign_mat, n_draws = 10000)

# Assign shock names to second dimension
dimnames(svar_sign_pre$irf_mean)[[2]]  <- colnames(sign_mat)
dimnames(svar_sign_pre$irf_lower)[[2]] <- colnames(sign_mat)
dimnames(svar_sign_pre$irf_upper)[[2]] <- colnames(sign_mat)

dimnames(svar_sign_post$irf_mean)[[2]]  <- colnames(sign_mat)
dimnames(svar_sign_post$irf_lower)[[2]] <- colnames(sign_mat)
dimnames(svar_sign_post$irf_upper)[[2]] <- colnames(sign_mat)

# --- 4d. Plot structural IRFs ---
plot_sign_irf <- function(svar_pre, svar_post,
                           impulse  = "NVDA_shock",
                           response = "SOXX",
                           n_ahead  = 5,
                           scale    = 100) {

  shock_idx    <- which(colnames(sign_mat) == impulse)
  response_idx <- which(rownames(sign_mat) == response)
  h <- 0:n_ahead

  pre_irf   <- svar_pre$irf_mean[response_idx,  shock_idx, ] * scale
  pre_lower <- svar_pre$irf_lower[response_idx, shock_idx, ] * scale
  pre_upper <- svar_pre$irf_upper[response_idx, shock_idx, ] * scale

  post_irf   <- svar_post$irf_mean[response_idx,  shock_idx, ] * scale
  post_lower <- svar_post$irf_lower[response_idx, shock_idx, ] * scale
  post_upper <- svar_post$irf_upper[response_idx, shock_idx, ] * scale

  ylim <- range(c(pre_lower, pre_upper, post_lower, post_upper))

  par(mfrow = c(1, 2))

  plot(h, pre_irf, type = "l", col = "steelblue", lwd = 2,
       ylim = ylim, xlab = "Weeks", ylab = "Response (%)",
       main = paste0("Sign SVAR: ", impulse, " \u2192 ", response, "\nPre-2015"))
  polygon(c(h, rev(h)), c(pre_upper, rev(pre_lower)),
          col = adjustcolor("steelblue", 0.2), border = NA)
  abline(h = 0, lty = 2, col = "gray50")

  plot(h, post_irf, type = "l", col = "coral", lwd = 2,
       ylim = ylim, xlab = "Weeks", ylab = "Response (%)",
       main = paste0("Sign SVAR: ", impulse, " \u2192 ", response, "\nPost-2015"))
  polygon(c(h, rev(h)), c(post_upper, rev(post_lower)),
          col = adjustcolor("coral", 0.2), border = NA)
  abline(h = 0, lty = 2, col = "gray50")

  par(mfrow = c(1, 1))
}

# Key shock-response pairs
plot_sign_irf(svar_sign_pre, svar_sign_post, "NVDA_shock",   "SOXX")
plot_sign_irf(svar_sign_pre, svar_sign_post, "Sector_shock", "NVDA")
plot_sign_irf(svar_sign_pre, svar_sign_post, "NVDA_shock",   "AMD")
plot_sign_irf(svar_sign_pre, svar_sign_post, "NVDA_shock",   "INTC")

# --- 4e. Forecast Error Variance Decomposition (FEVD) ---
# Note: IRF decays to near-zero by horizon 1 (fast-mean-reverting VAR(1)),
# so FEVD is effectively constant across horizons and reflects impact-period shares.

fevd_sign <- function(svar_result, sign_mat, n_ahead = 5) {

  n      <- dim(svar_result$irf_mean)[1]
  vnames <- dimnames(svar_result$irf_mean)[[1]]
  snames <- colnames(sign_mat)

  fevd_list <- list()

  for (resp_idx in seq_len(n)) {
    resp <- vnames[resp_idx]

    mse_total <- numeric(n_ahead + 1)
    mse_shock <- matrix(0, nrow = n_ahead + 1, ncol = n)

    for (h in 0:n_ahead) {
      for (s in seq_len(n)) {
        # Cumulative squared IRF from horizon 0 to h
        contrib <- sum(svar_result$irf_mean[resp_idx, s, 1:(h+1)]^2)
        mse_shock[h+1, s] <- contrib
      }
      mse_total[h+1] <- sum(mse_shock[h+1, ])
    }

    fevd_mat <- sweep(mse_shock, 1, mse_total, "/")
    colnames(fevd_mat) <- snames
    rownames(fevd_mat) <- 0:n_ahead

    fevd_list[[resp]] <- fevd_mat
  }

  return(fevd_list)
}

fevd_pre  <- fevd_sign(svar_sign_pre,  sign_mat, n_ahead = 5)
fevd_post <- fevd_sign(svar_sign_post, sign_mat, n_ahead = 5)

# Print summary at key horizons
cat("=== FEVD Pre-2015 ===\n")
for (v in c("NVDA", "AMD", "INTC", "SOXX")) {
  cat("\n", v, "(horizons 1, 3, 5):\n")
  print(round(fevd_pre[[v]][c(2, 4, 6), ], 3))
}

cat("\n=== FEVD Post-2015 ===\n")
for (v in c("NVDA", "AMD", "INTC", "SOXX")) {
  cat("\n", v, "(horizons 1, 3, 5):\n")
  print(round(fevd_post[[v]][c(2, 4, 6), ], 3))
}

# Impact-period summary (horizon 0)
cat("\n=== FEVD Summary (impact, horizon 0) ===\n")
for (v in c("NVDA", "AMD", "INTC", "SOXX")) {
  cat("\n", v, ":\n")
  df <- rbind(
    "Pre-2015"  = round(fevd_pre[[v]][1, ],  3),
    "Post-2015" = round(fevd_post[[v]][1, ], 3)
  )
  print(df)
}

# Plot FEVD comparison: pre vs post for each variable
plot_fevd_comparison <- function(fevd_pre, fevd_post, variable, horizon = 5) {

  pre  <- fevd_pre[[variable]][horizon + 1, ]
  post <- fevd_post[[variable]][horizon + 1, ]

  mat <- rbind(pre, post)
  rownames(mat) <- c("Pre-2015", "Post-2015")

  barplot(mat,
          beside    = TRUE,
          col       = c("steelblue", "coral"),
          main      = paste("FEVD:", variable, "| horizon", horizon),
          ylab      = "Share of variance",
          legend    = rownames(mat),
          las       = 2,
          cex.names = 0.85,
          ylim      = c(0, 1))
}

par(mfrow = c(2, 2))
plot_fevd_comparison(fevd_pre, fevd_post, "NVDA")
plot_fevd_comparison(fevd_pre, fevd_post, "AMD")
plot_fevd_comparison(fevd_pre, fevd_post, "INTC")
plot_fevd_comparison(fevd_pre, fevd_post, "SOXX")
par(mfrow = c(1, 1))

# ==============================================================================
# BLOCK 5: DCC-GJR-GARCH
# ==============================================================================

# --- Step 1: ARCH-LM test on VAR residuals ---
cat("\n===== ARCH-LM test on VAR residuals (pre-2015) =====\n")
resid_pre <- residuals(var_pre)
for (col in colnames(resid_pre)) {
  test <- ArchTest(resid_pre[, col], lags = 5)
  cat(sprintf("  %s: Chi2=%.3f, p-value=%.4f %s\n",
              col, test$statistic, test$p.value,
              ifelse(test$p.value < 0.05, "*** ARCH effects present", "")))
}

cat("\n===== ARCH-LM test on VAR residuals (post-2015) =====\n")
resid_post <- residuals(var_post)
for (col in colnames(resid_post)) {
  test <- ArchTest(resid_post[, col], lags = 5)
  cat(sprintf("  %s: Chi2=%.3f, p-value=%.4f %s\n",
              col, test$statistic, test$p.value,
              ifelse(test$p.value < 0.05, "*** ARCH effects present", "")))
}

# --- Step 2: Univariate GARCH specifications ---

# NVDA: GJR-GARCH to capture leverage effect
spec_nvda_bef <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "sstd"
)
spec_nvda_aft <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "sstd"
)

# AMD: standard sGARCH
spec_amd_bef <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "sstd"
)
spec_amd_aft <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "sstd"
)

# INTC: same spec for both periods
spec_int <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "sstd"
)

# SOXX
spec_soxx_bef <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "sstd"
)
spec_soxx_aft <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "sstd"
)

# --- Step 3: DCC specification ---
dcc_spec_bef <- dccspec(
  uspec        = multispec(list(spec_nvda_bef, spec_amd_bef, spec_int, spec_soxx_bef)),
  dccOrder     = c(1, 1),
  distribution = "mvt"
)

dcc_spec_aft <- dccspec(
  uspec        = multispec(list(spec_nvda_aft, spec_amd_aft, spec_int, spec_soxx_aft)),
  dccOrder     = c(1, 1),
  distribution = "mvt"
)

# --- Step 4: Fit DCC-GARCH ---
dcc_fit_pre  <- dccfit(dcc_spec_bef, data = df_pre,  fit.control = list(eval.se = TRUE))
dcc_fit_post <- dccfit(dcc_spec_aft, data = df_post, fit.control = list(eval.se = TRUE))

dcc_fit_pre
dcc_fit_post

# --- Step 5: Extract dynamic correlations ---
pairs_names <- c("NVDA-AMD", "NVDA-INTC", "NVDA-SOXX",
                 "AMD-INTC",  "AMD-SOXX",  "INTC-SOXX")

extract_dcc_corr <- function(dcc_fit, period_label) {
  R     <- rcor(dcc_fit)
  Tobs  <- dim(R)[3]
  pairs_idx <- list(
    c(2, 1), c(3, 1), c(4, 1),
    c(3, 2), c(4, 2), c(4, 3)
  )
  df_list <- lapply(seq_along(pairs_idx), function(k) {
    i <- pairs_idx[[k]][1]
    j <- pairs_idx[[k]][2]
    data.frame(
      t      = 1:Tobs,
      corr   = R[i, j, ],
      pair   = pairs_names[k],
      period = period_label
    )
  })
  do.call(rbind, df_list)
}

corr_pre  <- extract_dcc_corr(dcc_fit_pre,  "Before 2015")
corr_post <- extract_dcc_corr(dcc_fit_post, "After 2015")

# --- Step 6: DCC parameter comparison (alpha, beta, persistence) ---
get_dcc_pars <- function(fit_obj) {
  all_coefs <- coef(fit_obj)
  dcc_a <- all_coefs[grep("dcca1", names(all_coefs))]
  dcc_b <- all_coefs[grep("dccb1", names(all_coefs))]
  return(c(alpha = as.numeric(dcc_a), beta = as.numeric(dcc_b)))
}

compare_pars <- data.frame(
  Parameter = c("DCC Alpha (a)", "DCC Beta (b)"),
  Pre_SB    = get_dcc_pars(dcc_fit_pre),
  Post_SB   = get_dcc_pars(dcc_fit_post)
)
compare_pars <- rbind(compare_pars,
                      c("Persistence (a+b)",
                        sum(compare_pars$Pre_SB),
                        sum(compare_pars$Post_SB)))
print(compare_pars)

# --- Step 7: Summary table of mean correlations ---
tbl_pre <- corr_pre %>%
  group_by(pair) %>%
  summarise(
    mean_corr = round(mean(corr), 4),
    min_corr  = round(min(corr),  4),
    max_corr  = round(max(corr),  4)
  )
cat("\n===== Mean dynamic correlations (pre-2015) =====\n")
print(tbl_pre)

tbl_post <- corr_post %>%
  group_by(pair) %>%
  summarise(
    mean_corr = round(mean(corr), 4),
    min_corr  = round(min(corr),  4),
    max_corr  = round(max(corr),  4)
  )
cat("\n===== Mean dynamic correlations (post-2015) =====\n")
print(tbl_post)

cat("\n===== Change in mean correlation (post minus pre) =====\n")
tbl_diff        <- merge(tbl_pre, tbl_post, by = "pair", suffixes = c("_pre", "_post"))
tbl_diff$delta  <- round(tbl_diff$mean_corr_post - tbl_diff$mean_corr_pre, 4)
print(tbl_diff[, c("pair", "mean_corr_pre", "mean_corr_post", "delta")])

# --- Step 8: Plots ---

# Normalize time within each period for comparability
corr_all <- rbind(
  corr_pre  %>% group_by(pair, period) %>% mutate(t_norm = (t - min(t)) / (max(t) - min(t))),
  corr_post %>% group_by(pair, period) %>% mutate(t_norm = (t - min(t)) / (max(t) - min(t)))
)

# 8a. All pairs, both periods (faceted)
ggplot(corr_all, aes(x = t_norm, y = corr, color = period)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ pair, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Before 2015" = "#2196F3", "After 2015" = "#E53935")) +
  labs(
    title    = "DCC-GARCH: Dynamic correlations by pair",
    subtitle = "Normalized time (0 = period start, 1 = period end)",
    x        = "Normalized time",
    y        = "Conditional correlation",
    color    = "Period"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

# 8b. Focus: pairs involving NVDA
corr_nvda <- corr_all %>% filter(grepl("NVDA", pair))

ggplot(corr_nvda, aes(x = t_norm, y = corr, color = period)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ pair, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Before 2015" = "#2196F3", "After 2015" = "#E53935")) +
  labs(
    title = "DCC-GARCH: Dynamic correlations with NVDA",
    x     = "Normalized time",
    y     = "Conditional correlation",
    color = "Period"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

# 8c. Boxplot: distribution of correlations pre vs post
ggplot(corr_all, aes(x = pair, y = corr, fill = period)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Before 2015" = "#90CAF9", "After 2015" = "#EF9A9A")) +
  labs(
    title = "Distribution of conditional correlations: before vs after 2015",
    x     = NULL,
    y     = "Conditional correlation",
    fill  = "Period"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

# --- Step 9: DCC-GARCH parameter summary ---
cat("\n===== DCC-GARCH parameters (pre-2015) =====\n")
print(dcc_fit_pre)

cat("\n===== DCC-GARCH parameters (post-2015) =====\n")
print(dcc_fit_post)

# ==============================================================================
# BLOCK 6: ROLLING CORRELATION (52-week window)
# ==============================================================================

nvda_full_xts <- rbind(nvd_pre, nvd_post)
amd_full_xts  <- rbind(amd_pre, amd_post)
intc_full_xts <- rbind(int_pre, int_post)
soxx_full_xts <- rbind(sox_pre, sox_post)

roll_nvda_amd  <- rollapply(merge(nvda_full_xts, amd_full_xts),  width = 52, FUN = function(x) cor(x[, 1], x[, 2]), by.column = FALSE, align = "right")
roll_nvda_intc <- rollapply(merge(nvda_full_xts, intc_full_xts), width = 52, FUN = function(x) cor(x[, 1], x[, 2]), by.column = FALSE, align = "right")
roll_nvda_soxx <- rollapply(merge(nvda_full_xts, soxx_full_xts), width = 52, FUN = function(x) cor(x[, 1], x[, 2]), by.column = FALSE, align = "right")

df_roll <- data.frame(
  date      = index(roll_nvda_amd),
  NVDA_AMD  = as.numeric(roll_nvda_amd),
  NVDA_INTC = as.numeric(roll_nvda_intc),
  NVDA_SOXX = as.numeric(roll_nvda_soxx)
)

df_long <- pivot_longer(df_roll, cols = -date, names_to = "pair", values_to = "correlation")

ggplot(df_long, aes(x = date, y = correlation, color = pair)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2015-07-27"), linetype = "dashed", color = "black", linewidth = 1) +
  annotate("text", x = as.Date("2015-07-27"), y = 0.9, label = "Structural break (2015)", hjust = -0.1) +
  labs(
    title = "Rolling correlation with NVDA (52-week window)",
    x     = "Date",
    y     = "Correlation",
    color = "Pair"
  ) +
  theme_minimal()

# ==============================================================================
# BLOCK 7: DCC-GARCH DIAGNOSTICS
# ==============================================================================
# Evaluates whether the fitted DCC-GARCH adequately captures
# ARCH effects and autocorrelation in each series.

series_names <- c("NVDA", "AMD", "INTC", "SOXX")

run_diagnostics <- function(dcc_fit, period_label) {
  cat("\n", strrep("=", 60), "\n")
  cat(" Diagnostics:", period_label, "\n")
  cat(strrep("=", 60), "\n")

  # Standardized residuals: z = epsilon / sigma
  z <- residuals(dcc_fit) / sigma(dcc_fit)
  colnames(z) <- series_names

  for (s in series_names) {
    cat(sprintf("\n--- %s ---\n", s))

    # ARCH-LM: should be insignificant if GARCH captured volatility clustering
    arch <- ArchTest(z[, s], lags = 5)
    cat(sprintf("  ARCH-LM(5):      Chi2=%.3f, p=%.4f  %s\n",
                arch$statistic, arch$p.value,
                ifelse(arch$p.value < 0.05, "** ARCH remains", "OK")))

    # Ljung-Box on residuals: checks for remaining autocorrelation in mean
    lb <- Box.test(z[, s], lag = 10, type = "Ljung-Box")
    cat(sprintf("  Ljung-Box(10):   Chi2=%.3f, p=%.4f  %s\n",
                lb$statistic, lb$p.value,
                ifelse(lb$p.value < 0.05, "** autocorrelation remains", "OK")))

    # Ljung-Box on squared residuals: checks for remaining volatility clustering
    lb2 <- Box.test(z[, s]^2, lag = 10, type = "Ljung-Box")
    cat(sprintf("  Ljung-Box^2(10): Chi2=%.3f, p=%.4f  %s\n",
                lb2$statistic, lb2$p.value,
                ifelse(lb2$p.value < 0.05, "** volatility clustering remains", "OK")))

    # Jarque-Bera: checks if standardized residuals are normally distributed
    jb <- jarque.bera.test(z[, s])
    cat(sprintf("  Jarque-Bera:     Chi2=%.3f, p=%.4f  %s\n",
                jb$statistic, jb$p.value,
                ifelse(jb$p.value < 0.05, "** non-normal residuals", "OK")))
  }

  # DCC cross-residual check: products zi*zj should show no autocorrelation
  cat("\n--- DCC cross-residual Ljung-Box (zi * zj, lag=10) ---\n")
  pairs_diag <- list(
    c("NVDA", "AMD"), c("NVDA", "INTC"), c("NVDA", "SOXX"),
    c("AMD",  "INTC"), c("AMD",  "SOXX"), c("INTC", "SOXX")
  )
  for (p in pairs_diag) {
    lb_cross <- Box.test(z[, p[1]] * z[, p[2]], lag = 10, type = "Ljung-Box")
    cat(sprintf("  %s x %s: Chi2=%.3f, p=%.4f  %s\n",
                p[1], p[2],
                lb_cross$statistic, lb_cross$p.value,
                ifelse(lb_cross$p.value < 0.05, "** DCC misspecification?", "OK")))
  }
}

run_diagnostics(dcc_fit_pre,  "Pre-2015")
run_diagnostics(dcc_fit_post, "Post-2015")

# ==============================================================================
# BLOCK 8: DCC-GARCH FORECAST (12 weeks ahead)
# ==============================================================================
# Forecasts conditional volatility (sigma) and correlations (rho)
# based on the post-2015 model — the current regime.

n_ahead <- 12

fc_post <- dccforecast(dcc_fit_post, n.ahead = n_ahead)

# --- Volatility forecast ---
sigma_fc <- sigma(fc_post)   # matrix [n_ahead x 4]
colnames(sigma_fc) <- series_names

cat("\n===== Volatility forecast: post-2015 model (weekly sigma) =====\n")
print(round(sigma_fc, 6))

# Annualize: weekly sigma * sqrt(52)
sigma_ann <- sigma_fc * sqrt(52)
cat("\n===== Annualized volatility forecast (sigma * sqrt(52)) =====\n")
print(round(sigma_ann, 4))

# Plot forecasted weekly sigma
sigma_df <- as.data.frame(sigma_fc) %>%
  setNames(series_names) %>%
  mutate(horizon = 1:n_ahead) %>%
  pivot_longer(-horizon, names_to = "series", values_to = "sigma")

ggplot(sigma_df, aes(x = horizon, y = sigma, color = series)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 1:n_ahead) +
  labs(
    title    = "DCC-GARCH forecast: conditional volatility (post-2015 model)",
    subtitle = "12-week ahead forecast of weekly sigma",
    x        = "Weeks ahead",
    y        = "Conditional sigma",
    color    = "Series"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# --- Correlation forecast ---
R_fc <- rcor(fc_post)[[1]]   # array [4 x 4 x n_ahead]

pairs_idx <- list(
  c(2, 1), c(3, 1), c(4, 1),
  c(3, 2), c(4, 2), c(4, 3)
)

corr_fc_df <- do.call(rbind, lapply(seq_along(pairs_idx), function(k) {
  i <- pairs_idx[[k]][1]
  j <- pairs_idx[[k]][2]
  data.frame(
    horizon = 1:n_ahead,
    corr    = R_fc[i, j, ],
    pair    = pairs_names[k]
  )
}))

cat("\n===== Correlation forecast: post-2015 model =====\n")
print(
  corr_fc_df %>%
    group_by(pair) %>%
    summarise(
      h1  = round(corr[horizon == 1],       4),
      h6  = round(corr[horizon == 6],       4),
      h12 = round(corr[horizon == n_ahead], 4)
    )
)

# Plot forecasted correlations
ggplot(corr_fc_df, aes(x = horizon, y = corr, color = pair)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = 1:n_ahead) +
  labs(
    title    = "DCC-GARCH forecast: conditional correlations (post-2015 model)",
    subtitle = "12-week ahead forecast",
    x        = "Weeks ahead",
    y        = "Forecasted conditional correlation",
    color    = "Pair"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
