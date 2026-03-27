# ==============================================================================
# Structural Break Analysis: NVIDIA (NVDA)
# GARCH-adjusted series, Sup-F / Ave-F / Exp-F, CUSUM, MOSUM, RE, ME,
# Bai-Perron, and GARCH forecasts with vs without structural break
# ==============================================================================

library("quantmod")
library("TSA")
library("forecast")
library("urca")
library("tseries")
library("lmtest")
library("rugarch")
library("fGarch")
library("strucchangeRcpp")

# ==============================================================================
# DATA
# ==============================================================================

getSymbols("NVDA", from = "2001-07-09", to = Sys.Date(), periodicity = "weekly")

nvd_cls <- ts(NVDA$NVDA.Close, start = c(2001, 7), frequency = 52)
nvd_cls <- na.omit(nvd_cls)
length(nvd_cls)
ts.plot(nvd_cls)

# Log-returns
stat_nvd <- ts(diff(log(nvd_cls), differences = 1), start = c(2001, 7), frequency = 52)

ts.plot(stat_nvd)

# Stationarity tests
adf.test(stat_nvd)
kpss.test(stat_nvd)
pp.test(stat_nvd)

# ==============================================================================
# ARMA + GARCH VOLATILITY ADJUSTMENT
# ==============================================================================

# ACF/PACF inspection — NVDA: AR(1) sufficient
acf(stat_nvd)
pacf(stat_nvd)
ar0model <- arima(stat_nvd, c(1, 0, 0))
coeftest(ar0model)

# Residual diagnostics
acf(ar0model$residuals)
Box.test(ar0model$residuals)
acf(ar0model$residuals^2)
Box.test(ar0model$residuals^2)   # ARCH effects present

# Fit GJR-GARCH(1,1) to capture leverage effect
spec <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE, archm = FALSE),
  distribution.model = "sstd"
)
garch_fit <- ugarchfit(spec, stat_nvd)
garch_fit

# GARCH-adjusted series (standardised residuals)
stat_nvd_adj <- stat_nvd / garch_fit@fit$sigma
ts.plot(stat_nvd_adj,
        main = "NVDA log-returns adjusted for conditional volatility",
        ylab = "y_adj")

# ==============================================================================
# STRUCTURAL BREAK TESTS ON GARCH-ADJUSTED SERIES
# ==============================================================================

# --- 1. Sup-F / Ave-F / Exp-F (Andrews / Andrews-Ploberger) ---

# Sup-F
stat_supf <- Fstats(stat_nvd_adj ~ 1, from = 0.1)
plot(stat_supf, alpha = 0.01, main = "Sup-F test — NVDA AR(0)")
lines(breakpoints(stat_supf))
breakpoints(stat_supf)
sctest(stat_supf, type = "supF")

# Ave-F
stat_avef <- Fstats(stat_nvd_adj ~ 1, from = 0.2)
plot(stat_avef, alpha = 0.01, aveF = TRUE, main = "Ave-F test — NVDA AR(0)")
lines(breakpoints(stat_avef))
breakpoints(stat_avef)
sctest(stat_avef, type = "aveF")

# Exp-F
sctest(stat_avef, type = "expF")

# --- 2. CUSUM tests ---

# OLS-CUSUM
stat_cusum_ols <- efp(stat_nvd_adj ~ 1, type = "OLS-CUSUM")
plot(stat_cusum_ols, alpha = 0.1, functional = NULL, main = "OLS-CUSUM — NVDA AR(0)")
sctest(stat_cusum_ols)

# Recursive CUSUM
stat_cusum_rec <- efp(stat_nvd_adj ~ 1, type = "Rec-CUSUM")
plot(stat_cusum_rec, alpha = 0.1, main = "Rec-CUSUM — NVDA AR(0)")
sctest(stat_cusum_rec)

# --- 3. MOSUM tests ---

# OLS-MOSUM
stat_mosum_ols <- efp(stat_nvd_adj ~ 1, h = 0.5, type = "OLS-MOSUM")
plot(stat_mosum_ols, alpha = 0.1, main = "OLS-MOSUM — NVDA AR(0)")
sctest(stat_mosum_ols)

# Recursive MOSUM
stat_mosum_rec <- efp(stat_nvd_adj ~ 1, h = 0.5, type = "Rec-MOSUM")
plot(stat_mosum_rec, alpha = 0.1, main = "Rec-MOSUM — NVDA AR(0)")
sctest(stat_mosum_rec)

# --- 4. Parameter stability: RE and ME ---

# RE (Recursive Estimates)
stat_re <- efp(stat_nvd_adj ~ 1, type = "RE")
plot(stat_re, alpha = 0.1, functional = NULL, main = "RE test — NVDA AR(0)")
sctest(stat_re)

# ME (Moving Estimates)
stat_me <- efp(stat_nvd_adj ~ 1, h = 0.5, type = "ME")
plot(stat_me, alpha = 0.1, functional = NULL, main = "ME test — NVDA AR(0)")
sctest(stat_me)

# --- 5. Partial structural break (intercept only) ---
stat_partial_const <- gefp(stat_nvd_adj ~ 1, parm = 1)
plot(stat_partial_const, alpha = 0.1, main = "Partial SB: intercept — NVDA ARMA(0,0)")
sctest(stat_partial_const)

# --- 6. Bai-Perron multiple breakpoint test ---

cat("\n===== Bai-Perron: raw log-returns =====\n")
bp_bai_orig <- breakpoints(stat_nvd ~ 1)
summary(bp_bai_orig)
plot(bp_bai_orig, main = "Bai-Perron — NVDA raw returns")
breakdates(bp_bai_orig)

cat("\n===== Bai-Perron: GARCH-adjusted series =====\n")
bp_bai_adj <- breakpoints(stat_nvd_adj ~ 1)
summary(bp_bai_adj)
plot(bp_bai_adj, main = "Bai-Perron — NVDA GARCH-adjusted")
breakdates(bp_bai_adj)

# ==============================================================================
# GARCH FORECASTS: FULL SAMPLE vs POST-BREAK SUBSAMPLE
# Structural break: week index 811 (July 2015)
# ==============================================================================

h        <- 52    # forecast horizon (weeks)
sb_index <- 811   # break index in stat_nvd

# --- Model 1: GJR-GARCH on full sample (ignoring break) ---
spec_full <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE, archm = FALSE),
  distribution.model = "sstd"
)
garch_full     <- ugarchfit(spec_full, stat_nvd[1:sb_index])
forecast_full  <- ugarchforecast(garch_full, n.ahead = h)
mu_full        <- as.numeric(fitted(forecast_full))
sigma_full     <- as.numeric(sigma(forecast_full))

# --- Model 2: GJR-GARCH on post-break subsample ---
stat_nvd_after <- stat_nvd[(sb_index + 1):length(stat_nvd)]
cat("Post-break subsample length:", length(stat_nvd_after), "weeks\n")

spec_sb    <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE, archm = FALSE),
  distribution.model = "sstd"
)
garch_sb    <- ugarchfit(spec_sb, stat_nvd_after)
forecast_sb <- ugarchforecast(garch_sb, n.ahead = h)
mu_sb       <- as.numeric(fitted(forecast_sb))
sigma_sb    <- as.numeric(sigma(forecast_sb))

# --- Parameter comparison ---
cat("\n===== Model 1: full sample =====\n")
cat("mu:", coef(garch_full)["mu"], "\n")
cat("1-step mean forecast:", mu_full[1], "\n")
cat("1-step volatility forecast:", sigma_full[1], "\n")

cat("\n===== Model 2: post-break (2015) =====\n")
cat("mu:", coef(garch_sb)["mu"], "\n")
cat("1-step mean forecast:", mu_sb[1], "\n")
cat("1-step volatility forecast:", sigma_sb[1], "\n")

# --- Price forecasts ---
last_price <- as.numeric(tail(NVDA$NVDA.Close, 1))
cat("Last observed price: $", last_price, "\n")

price_full <- last_price * exp(cumsum(mu_full))
price_sb   <- last_price * exp(cumsum(mu_sb))

# ==============================================================================
# PLOTS
# ==============================================================================

n_show   <- 100
last_obs <- as.numeric(tail(na.omit(nvd_cls), n_show))
time_obs <- (-n_show + 1):0
time_fct <- 1:h

# --- Plot 1: mean return forecast ---
par(mfrow = c(1, 1))
plot(1:h, mu_full, type = "l", col = "blue", lwd = 2,
     ylim = range(c(mu_full, mu_sb)),
     xlab = "Weeks ahead", ylab = "Log-return",
     main = "NVDA: mean return forecast\nFull sample (blue) vs post-break 2015 (red)")
lines(1:h, mu_sb, col = "red", lwd = 2)
abline(h = 0, lty = 2, col = "grey")
legend("topright",
       legend = c("No break (full sample)", "With break (post-2015)"),
       col = c("blue", "red"), lwd = 2)

# --- Plot 2: volatility forecast ---
plot(1:h, sigma_full, type = "l", col = "blue", lwd = 2,
     ylim = range(c(sigma_full, sigma_sb)),
     xlab = "Weeks ahead", ylab = "Conditional std. deviation",
     main = "NVDA: volatility forecast\nFull sample (blue) vs post-break 2015 (red)")
lines(1:h, sigma_sb, col = "red", lwd = 2)
legend("topright",
       legend = c("No break (full sample)", "With break (post-2015)"),
       col = c("blue", "red"), lwd = 2)

# --- Plot 3: price forecast ---
plot(time_obs, last_obs, type = "l", col = "black", lwd = 1.5,
     xlim = c(-n_show, h), ylim = range(c(0, 350)),
     xlab = "Weeks (0 = last observation)", ylab = "NVDA stock price ($)",
     main = "NVDA: 52-week price forecast\nFull sample (blue) vs post-break 2015 (red)")
lines(time_fct, price_full, col = "blue", lwd = 2)
lines(time_fct, price_sb,   col = "red",  lwd = 2)
abline(v = 0, lty = 2, col = "grey")
legend("topleft",
       legend = c("Observed", "No break (full sample)", "With break (post-2015)"),
       col = c("black", "blue", "red"), lwd = c(1.5, 2, 2))

# ==============================================================================
# SIMULATED FORECAST PATHS
# ==============================================================================

n_sim <- 15

sim_full <- ugarchsim(garch_full, n.sim = h, m.sim = n_sim, rseed = 42)
sim_sb   <- ugarchsim(garch_sb,   n.sim = h, m.sim = n_sim, rseed = 42)

paths_full <- fitted(sim_full)
paths_sb   <- fitted(sim_sb)

price_paths_full <- matrix(NA, nrow = h, ncol = n_sim)
price_paths_sb   <- matrix(NA, nrow = h, ncol = n_sim)
for (i in 1:n_sim) {
  price_paths_full[, i] <- last_price * exp(cumsum(paths_full[, i]))
  price_paths_sb[, i]   <- last_price * exp(cumsum(paths_sb[, i]))
}

# Side-by-side: price paths
par(mfrow = c(1, 2))

plot(time_obs, last_obs, type = "l", col = "black", lwd = 2,
     xlim = c(-n_show, h),
     ylim = range(c(last_obs, price_paths_full, price_full)),
     xlab = "Weeks", ylab = "Price ($)",
     main = "Model 1: full sample\n15 simulated paths")
for (i in 1:n_sim) lines(time_fct, price_paths_full[, i], col = adjustcolor("blue", 0.3), lwd = 1)
lines(time_fct, price_full, col = "blue", lwd = 2.5)
abline(v = 0, lty = 2, col = "grey")
legend("topleft",
       legend = c("Observed", "Mean forecast", "Simulated paths"),
       col = c("black", "blue", adjustcolor("blue", 0.3)), lwd = c(2, 2.5, 1))

plot(time_obs, last_obs, type = "l", col = "black", lwd = 2,
     xlim = c(-n_show, h),
     ylim = range(c(last_obs, price_paths_sb, price_sb)),
     xlab = "Weeks", ylab = "Price ($)",
     main = "Model 2: post-break 2015\n15 simulated paths")
for (i in 1:n_sim) lines(time_fct, price_paths_sb[, i], col = adjustcolor("red", 0.3), lwd = 1)
lines(time_fct, price_sb, col = "red", lwd = 2.5)
abline(v = 0, lty = 2, col = "grey")
legend("topleft",
       legend = c("Observed", "Mean forecast", "Simulated paths"),
       col = c("black", "red", adjustcolor("red", 0.3)), lwd = c(2, 2.5, 1))

par(mfrow = c(1, 1))

# Combined price paths
plot(time_obs, last_obs, type = "l", col = "black", lwd = 2,
     xlim = c(-n_show, h), ylim = range(c(last_obs, 370)),
     xlab = "Weeks (0 = last observation)", ylab = "NVDA stock price ($)",
     main = "NVDA simulated price paths:\nFull sample (blue) vs post-break 2015 (red)")
for (i in 1:n_sim) {
  lines(time_fct, price_paths_full[, i], col = adjustcolor("blue", 0.25), lwd = 1)
  lines(time_fct, price_paths_sb[, i],   col = adjustcolor("red",  0.25), lwd = 1)
}
lines(time_fct, price_full, col = "blue", lwd = 2.5)
lines(time_fct, price_sb,   col = "red",  lwd = 2.5)
abline(v = 0, lty = 2, col = "grey")
legend("topleft",
       legend = c("Observed", "Mean forecast (full)", "Mean forecast (post-break)",
                  "Paths (full)", "Paths (post-break)"),
       col = c("black", "blue", "red", adjustcolor("blue", 0.4), adjustcolor("red", 0.4)),
       lwd = c(2, 2.5, 2.5, 1, 1))

# ==============================================================================
# SIMULATED VOLATILITY PATHS
# ==============================================================================

sigma_paths_full <- sigma(sim_full)
sigma_paths_sb   <- sigma(sim_sb)

par(mfrow = c(1, 2))

plot(1:h, sigma_paths_full[, 1], type = "l",
     col = adjustcolor("blue", 0.3), lwd = 1,
     ylim = range(c(sigma_paths_full, sigma_full)),
     xlab = "Weeks ahead", ylab = "Conditional std. deviation",
     main = "Volatility: full sample\n15 simulated paths")
for (i in 2:n_sim) lines(1:h, sigma_paths_full[, i], col = adjustcolor("blue", 0.3), lwd = 1)
lines(1:h, sigma_full, col = "blue", lwd = 2.5)
legend("topright",
       legend = c("Mean forecast", "Simulated paths"),
       col = c("blue", adjustcolor("blue", 0.3)), lwd = c(2.5, 1))

plot(1:h, sigma_paths_sb[, 1], type = "l",
     col = adjustcolor("red", 0.3), lwd = 1,
     ylim = range(c(sigma_paths_sb, sigma_sb)),
     xlab = "Weeks ahead", ylab = "Conditional std. deviation",
     main = "Volatility: post-break 2015\n15 simulated paths")
for (i in 2:n_sim) lines(1:h, sigma_paths_sb[, i], col = adjustcolor("red", 0.3), lwd = 1)
lines(1:h, sigma_sb, col = "red", lwd = 2.5)
legend("topright",
       legend = c("Mean forecast", "Simulated paths"),
       col = c("red", adjustcolor("red", 0.3)), lwd = c(2.5, 1))

par(mfrow = c(1, 1))

# Combined volatility paths
plot(1:h, sigma_full, type = "l", col = "blue", lwd = 2.5,
     ylim = c(0.05, 0.08),
     xlab = "Weeks ahead", ylab = "Conditional std. deviation",
     main = "NVDA volatility forecast:\nFull sample (blue) vs post-break 2015 (red)")
for (i in 1:n_sim) {
  lines(1:h, sigma_paths_full[, i], col = adjustcolor("blue", 0.25), lwd = 1)
  lines(1:h, sigma_paths_sb[, i],   col = adjustcolor("red",  0.25), lwd = 1)
}
lines(1:h, sigma_full, col = "blue", lwd = 2.5)
lines(1:h, sigma_sb,   col = "red",  lwd = 2.5)
legend("topright",
       legend = c("Mean forecast (full)", "Mean forecast (post-break)",
                  "Paths (full)", "Paths (post-break)"),
       col = c("blue", "red", adjustcolor("blue", 0.4), adjustcolor("red", 0.4)),
       lwd = c(2.5, 2.5, 1, 1))
