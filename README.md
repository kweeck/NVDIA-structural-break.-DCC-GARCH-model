# NVIDIA Structural Break & Semiconductor Sector Dynamics

![Language](https://img.shields.io/badge/language-R-276DC3?style=flat-square&logo=r)
![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)
![Status](https://img.shields.io/badge/status-complete-brightgreen?style=flat-square)
![SSRN](https://img.shields.io/badge/SSRN-working%20paper-orange?style=flat-square)

> Independent research project — applied financial econometrics using structural break detection, VAR, Granger causality, IRF, and DCC-GJR-GARCH on 25 years of semiconductor sector data.

---

## Research Question

Did NVIDIA's 2015 structural break alter not only its return level, but also its **role within the system of interdependencies** in the semiconductor sector?

---

## Key Findings

| Finding | Before 2015 | After 2015 |
|---|---|---|
| Granger causality NVDA → SOXX | ✅ significant (p = 0.026) | ❌ disappeared (p = 0.219) |
| Granger causality SOXX → NVDA | ❌ absent | ✅ emerged (p = 0.010) |
| NVDA momentum (own lag) | ❌ absent | ✅ significant |
| NVDA mean return drift | near zero | positive & significant |
| DCC reaction speed (dcca1) | 0.013 | 0.034 |
| DCC persistence (dccb1) | 0.972 | 0.908 |
| NVDA–AMD correlation | ~0.50 | ~0.60 ↑ |
| NVDA–INTC correlation | ~0.50 | ~0.40 ↓ |

**Summary:** NVIDIA transitioned from a niche sector leader whose lagged returns predicted the broad index, to an autonomous AI-driven asset subject to institutional rotation dynamics. The sector split into two camps — NVDA/AMD converged as AI companies, while INTC diverged.

---

## Repository Structure

```
├── R/
│   ├── structural_break/
│   │   ├── sb_nvda.R       # Bai-Perron & CUSUM for NVDA
│   │   ├── sb_amd.R        # Structural break analysis for AMD
│   │   ├── sb_intc.R       # Structural break analysis for INTC
│   │   └── sb_soxx.R       # Structural break analysis for SOXX
│   └── analysis/
│       └── var_dcc.R       # VAR, Granger, IRF, DCC-GARCH, rolling correlation
├── images/                 # Plots referenced in this README
└── report/
    └── report.Rmd          # Full report (Rmd)
```

---

## Data

| Parameter | Value |
|---|---|
| Source | Yahoo Finance via `quantmod` |
| Frequency | Weekly log-returns |
| Period | 2001 – present |
| Assets | NVDA, AMD, INTC, SOXX ETF |
| Breakpoint | July 2015 |
| Pre-break obs. | 733 |
| Post-break obs. | 555 |

---

## Methodology

### Stage 1 — Structural Break Detection

Applied to GARCH-adjusted returns (AR(1) + GJR-GARCH(1,1) with skewed-t innovations) to remove volatility clustering before testing for mean shifts.

| Test | Result |
|---|---|
| Ave-F / Exp-F | Significant |
| OLS-CUSUM | Significant break in 2015 |
| Recursive Estimates | Significant break in 2015 |
| Partial structural break | Significant break in 2015 |
| Bai-Perron |Insignificant, indicates to 2015 |

![Ave-F](images/ave.png)
![OLS-CUSUM](images/olscums.png)
![RE](images/re.png)
![Partial](images/part.png)

---

### Stage 2 — Forecasting Impact

Comparing GARCH models estimated on the full sample vs. the post-break subsample shows materially different return and volatility forecasts — ignoring the break leads to systematic underestimation of post-2015 drift and volatility persistence.

![Price Forecast](images/pricefore.png)
![Volatility Forecast](images/volatil.png)

---

### Stage 3 — VAR Model

VAR(1) selected unanimously by AIC, HQ, SC, FPE for both sub-periods. Consistent with weak-form market efficiency at weekly frequency.

Residual diagnostics (Hosking, Li-McLeod) confirm white noise residuals. CUSUM stability confirms the 2015 split is methodologically justified.

**Pre-2015:** Only one significant coefficient in the entire VAR matrix — NVDA.l1 in the SOXX equation (β = 0.064, p = 0.029). NVDA led the sector through lags.

**Post-2015:** NVDA acquires a significant own lag (β = 0.158, p = 0.020) and positive drift (β = 0.010, p = 0.0002). The lagged leadership disappears.

---

### Stage 4 — Granger Causality

**System-level tests:**

| Variable | Pre-2015 (p) | Post-2015 (p) |
|---|---|---|
| NVDA → others | **0.037** ✅ | 0.411 ❌ |
| AMD → others | 0.580 | 0.921 |
| INTC → others | 0.161 | 0.509 |
| SOXX → others | 0.302 | 0.306 |

**Pairwise tests (key directions):**

| Direction | Pre-2015 (p) | Post-2015 (p) |
|---|---|---|
| NVDA → SOXX | **0.026** | 0.219 |
| SOXX → NVDA | 0.953 | **0.010** |
| NVDA → AMD | 0.103 | 0.070 |

The reversal of SOXX → NVDA causality is interpreted as **institutional rotation**: large funds hold NVDA as a semiconductor bet and rebalance relative to the ETF benchmark.

---

### Stage 5 — Impulse Response Functions

**NVDA → SOXX:** Significant positive response at week 1 before 2015 (~0.065, CI excludes zero). After 2015: response shrinks (~0.047) and becomes statistically insignificant.

**SOXX → NVDA:** No significant response before 2015. After 2015: persistently negative direction (~−0.2), consistent with Granger p = 0.010, interpreted as rotation effect.

![IRF: NVDA → SOXX](images/irfnvd.png)
![IRF: SOXX → NVDA](images/irfsoxx.png)

---

### Stage 6 — DCC-GJR-GARCH

Univariate GJR-GARCH models (skewed-t innovations) estimated per asset, then combined in a DCC(1,1) structure. Estimated separately for each sub-period.

**NVDA post-2015:** beta1 = 0.996 — volatility process approaches a near-random walk. Large AI-driven moves dissipate very slowly.

**Dynamic correlation changes:**

| Pair | Pre-2015 | Post-2015 | Interpretation |
|---|---|---|---|
| NVDA–AMD | ~0.50 | ~0.60 ↑ | Converged as AI companies |
| NVDA–INTC | ~0.50 | ~0.40 ↓ | Diverged sharply |
| NVDA–SOXX | ~0.70 | ~0.75 ↑ | Deeper sector integration |
| AMD–INTC | ~0.50 | ~0.30 ↓ | AMD moved to AI/gaming |

![Dynamic Correlations](images/dcorrall.png)
![Boxplot](images/boxplot.png)

---

### Stage 7 — Rolling Correlation (Robustness Check)

52-week rolling window provides a non-parametric benchmark confirming DCC results. Post-2015 divergence of NVDA–INTC and relative stability of NVDA–AMD are clearly visible.

![Rolling Correlation](images/rollcor.png)

---

## Reproducibility

Data is downloaded automatically via `quantmod`. No local files required.

```r
install.packages(c(
  "quantmod", "TSA", "forecast", "urca", "tseries",
  "lmtest", "rugarch", "fGarch", "strucchange",
  "strucchangeRcpp", "vars", "rmgarch",
  "ggplot2", "dplyr", "tidyr", "zoo", "FinTS", "portes"
))
```

**Run order:**

```r
# Stage 1: Structural break analysis
source("R/structural_break/sb_nvda.R")

# Stage 2: Multivariate analysis
source("R/analysis/var_dcc.R")
```

## Limitations

- The breakpoint (July 2015) is treated as fixed based on test consensus
- Weekly frequency may smooth intra-week dynamics
- GARCH orders selected via standard diagnostics; exhaustive grid search not conducted
- DCC assumes a constant unconditional correlation matrix within each sub-period


## License

MIT — see [LICENSE](LICENSE) for details.
