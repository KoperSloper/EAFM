# Load required libraries
library(tidyverse)
library(dbscan)
library(PerformanceAnalytics)
library(xts)
library(covglasso)
library(HierPortfolios)
library(rugarch)    # for ugarchspec & related functions
library(rmgarch)    # for dccspec & dccroll

# -------------------------
# 1. Prepare and clean data
# -------------------------

# Get full date sequence and count
all_dates  <- sort(unique(even_more_stocks_realized$date))
n_dates    <- length(all_dates)

# Identify tickers (PERMNOs) with invalid returns or missing bid/ask
bad_permnos <- even_more_stocks_realized %>%
  filter(
    !str_detect(RET, "^[-]?[0-9]+(\\.[0-9]+)?$") |
      is.na(BIDLO) |
      is.na(ASKHI)
  ) %>%
  pull(PERMNO) %>%
  unique()

# Keep only those tickers that have no missing/bad data for every date
good_permnos <- even_more_stocks_realized %>%
  filter(!PERMNO %in% bad_permnos) %>%
  group_by(PERMNO) %>%
  summarize(n_obs = n_distinct(date), .groups = "drop") %>%
  filter(n_obs == n_dates) %>%
  pull(PERMNO)

# Subset to “good” tickers and compute realized volatility (Parkinson)
cleaned <- even_more_stocks_realized %>%
  filter(PERMNO %in% good_permnos) %>%
  mutate(
    RET      = as.numeric(RET),
    date     = as.Date(date),
    park_var = (log(ASKHI / BIDLO)^2) / (4 * log(2)),
    park_vol = sqrt(park_var)
  )

# Convert to wide matrices: returns and realized volatility
returns_wide <- cleaned %>%
  select(PERMNO, date, RET) %>%
  pivot_wider(names_from = PERMNO, values_from = RET) %>%
  arrange(date)

range_wide <- cleaned %>%
  select(PERMNO, date, park_vol) %>%
  pivot_wider(names_from = PERMNO, values_from = park_vol) %>%
  arrange(date)

# Extract date vector and numeric matrices
dates       <- returns_wide$date
returns_mat <- returns_wide %>% select(-date) %>% as.matrix()
range_mat   <- range_wide  %>% select(-date) %>% as.matrix()

rownames(returns_mat) <- as.character(dates)
rownames(range_mat)   <- as.character(dates)

p_assets <- ncol(returns_mat)
N        <- nrow(returns_mat)

# -----------------------------------
# 2. Fit a DCC‐GARCH(1,1) on returns
# -----------------------------------

# Specify Realized‐GARCH(1,1) univariate model
spec_rg <- ugarchspec(
  variance.model = list(model = "realGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

# Build a multispec for all assets
uspec_rg  <- multispec(replicate(p_assets, spec_rg))

dcc_spec_rg <- dccspec(
  uspec        = uspec_rg,
  dccOrder     = c(1, 1),
  model        = "aDCC",
  distribution = "mvnorm"
)

# Rolling‐window DCC estimation
rolling_rg <- dccroll(
  spec            = dcc_spec_rg,
  data            = returns_mat,
  realizedVol     = xts(range_mat, order.by = dates),
  n.ahead         = 1,
  forecast.length = 1000,
  refit.every     = 5,
  refit.window    = "moving",
  window.size     = 1000,
  fit.control     = list(eval.se = FALSE)
)

# Extract conditional standard deviations, covariances, correlations
sigma_mat   <- sigma(rolling_rg)          # T×p_assets
cov_array   <- rcov(rolling_rg)           # p_assets×p_assets×T
corr_array  <- rcor(rolling_rg)           # p_assets×p_assets×T

n_times      <- dim(cov_array)[3]
asset_names  <- colnames(sigma_mat)

# -------------------------------------------
# 3. Sparse covariance via graphical LASSO
# -------------------------------------------

eps <- 1e-6

sparse_cov_array  <- array(NA, dim = c(p_assets, p_assets, n_times))
sparse_corr_array <- array(NA, dim = c(p_assets, p_assets, n_times))

for (t in seq_len(n_times)) {
  S_t   <- cov_array[ , , t]
  S_reg <- S_t + diag(eps, p_assets)    # ensure positive definiteness
  fit   <- covglasso(S = S_reg, n = 1000, lambda = NULL, path = FALSE)
  sparse_cov_array[ , , t]  <- fit$sigma
  sparse_corr_array[ , , t] <- cov2cor(fit$sigma)
}

# --------------------------------------------------
# 4. “Naïve” sparse covariance (rolling empirical)
# --------------------------------------------------

naive_sparse_cov_array  <- array(NA, dim = c(p_assets, p_assets, n_times))
naive_sparse_corr_array <- naive_sparse_cov_array

for (t in seq_len(n_times)) {
  ret_row    <- (N - n_times) + t
  if (ret_row <= 1000) stop("Need ≥1000 prior obs for naive cov")
  window_idx <- (ret_row - 1000):(ret_row - 1)
  S_emp      <- cov(returns_mat[window_idx, ], use = "pairwise.complete.obs")
  S_reg      <- S_emp + diag(eps, p_assets)
  fit        <- covglasso(S = S_reg, n = length(window_idx), lambda = NULL, path = FALSE)
  naive_sparse_cov_array[,,t]  <- fit$sigma
  naive_sparse_corr_array[,,t] <- cov2cor(fit$sigma)
}

# ------------------------------------------------
# 5. Estimate spread per stock (H‐L spread method)
# ------------------------------------------------

# Define spread estimator function
HLSpreadEstimatorImposePositive <- function(highs, lows) {
  beta  <- (log(highs[1] / lows[1]))^2 + (log(highs[2] / lows[2]))^2
  H     <- max(highs);  L <- min(lows)
  gamma <- (log(H / L))^2
  alpha <- (sqrt(2 * beta) - sqrt(beta)) / (3 - 2 * sqrt(2)) -
    sqrt(gamma / (3 - 2 * sqrt(2)))
  s     <- 2 * (exp(alpha) - 1) / (1 + exp(alpha))
  max(s, 0)
}

# Compute cross‐sectional spread by stock
spread_estimates <- even_more_stocks_realized %>%
  filter(PERMNO %in% good_permnos) %>%
  arrange(PERMNO, date) %>%
  group_by(PERMNO) %>%
  mutate(
    lag_bid   = lag(BIDLO),
    lag_ask   = lag(ASKHI),
    spread_cs = if_else(
      !is.na(lag_bid),
      mapply(
        HLSpreadEstimatorImposePositive,
        c(ASKHI, lag_ask),
        c(BIDLO, lag_bid)
      ),
      NA_real_
    )
  ) %>%
  ungroup()

spread_summary_by_stock <- spread_estimates %>%
  filter(!is.na(spread_cs)) %>%
  group_by(PERMNO) %>%
  summarise(
    n_obs         = n(),
    mean_spread   = mean(spread_cs),
    median_spread = median(spread_cs),
    .groups = "drop"
  ) %>%
  mutate(
    c_spread = if_else(median_spread > 0, median_spread / 2, mean_spread / 2)
  )

spread_vec <- spread_summary_by_stock %>%
  arrange(match(PERMNO, good_permnos)) %>%
  pull(c_spread)

# ------------------------------------------------------
# 6. Build price and volume matrices, check consistency
# ------------------------------------------------------

filtered_pv <- even_more_stocks_realized_volume_price %>%
  filter(PERMNO %in% good_permnos)

price_mat_df <- filtered_pv %>%
  select(date, PERMNO, PRC) %>%
  pivot_wider(names_from = PERMNO, values_from = PRC) %>%
  arrange(date)

vol_mat_df <- filtered_pv %>%
  select(date, PERMNO, VOL) %>%
  pivot_wider(names_from = PERMNO, values_from = VOL) %>%
  arrange(date)

# Extract numeric matrices
dates_price <- price_mat_df$date
price_mat   <- price_mat_df %>% select(-date) %>% as.matrix()
vol_mat     <- vol_mat_df   %>% select(-date) %>% as.matrix()

# Consistency checks
stopifnot(
  nrow(price_mat) == nrow(returns_mat),
  ncol(price_mat) == ncol(returns_mat),
  all(colnames(price_mat) == colnames(returns_mat))
)

# -----------------------------------------------
# 7. Set up backtest parameters and containers
# -----------------------------------------------

validation_len <- 50
test_block     <- 10
tau_grid       <- seq(0.20, 0.35, by = 0.05)
threshold_turn <- 0.8    
c_i            <- 0.8
beta           <- 1.5
initial_aum    <- 1e6

# Preallocate tracking vectors
rets_dhrp   <- rets_naive <- rets_ew  <- numeric(n_times)
cost_dhrp   <- cost_naive <- cost_ew   <- numeric(n_times)
aum_dhrp    <- aum_naive  <- aum_ew    <- numeric(n_times)

# Initial weights and AUM for each strategy
w_dhrp_prev   <- rep(1 / p_assets, p_assets)
w_naive_prev  <- rep(1 / p_assets, p_assets)
w_ew_prev     <- rep(1 / p_assets, p_assets)

aum_dhrp[1:validation_len]  <- initial_aum
aum_naive[1:validation_len] <- initial_aum
aum_ew[1:validation_len]    <- initial_aum

# --------------------------------------------
# 8. Define helper functions for rebalancing
# --------------------------------------------

compute_w_dhrp <- function(Sigma_t, tau, t_index) {
  tryCatch({
    if (any(!is.finite(Sigma_t))) {
      return(rep(1 / p_assets, p_assets))
    }
    eig <- eigen(Sigma_t, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig <= 0)) {
      Sigma_t <- Sigma_t + diag(1e-5, p_assets)
    }
    DHRP_Portfolio(covar = Sigma_t, tau = tau,
                   UB = rep(1, p_assets), LB = rep(0, p_assets))$w
  }, error = function(e) {
    rep(1 / p_assets, p_assets)
  })
}

compute_w_naive <- function(Sigma_t) {
  w <- tryCatch(solve(Sigma_t, rep(1, p_assets)),
                error = function(e) rep(1, p_assets))
  w / sum(w)
}

calculate_transaction_costs <- function(w_new, w_prev, aum_prev,
                                        prices, volumes, sigma_diag,
                                        spread_vec) {
  Q            <- abs(w_new - w_prev) * aum_prev / prices
  spread_costs <- spread_vec * Q * prices
  q_today      <- volumes * prices
  d_tilde      <- Q * prices
  a_vec        <- (c_i * sigma_diag) / (q_today^(beta - 1))
  impact_costs <- a_vec * abs(d_tilde)^beta
  total_costs  <- spread_costs + impact_costs
  list(
    per_asset_cost = total_costs,
    portfolio_cost = sum(total_costs, na.rm = TRUE)
  )
}

# --------------------------------------------
# 9. Main backtest loops: validation + testing
# --------------------------------------------

day_counter <- 0

for (t0 in seq(validation_len + 1, n_times, by = test_block)) {
  # --- Validation to pick best tau ---
  val_idx  <- (t0 - validation_len):(t0 - 1)
  best_sh  <- -Inf
  best_tau <- NA_real_
  
  for (tau in tau_grid) {
    rets_val   <- numeric(length(val_idx))
    prev_d     <- rep(1 / p_assets, p_assets)
    prev_aum   <- if (t0 == validation_len + 1) initial_aum else aum_dhrp[t0 - 1]
    
    for (i in seq_along(val_idx)) {
      s      <- val_idx[i]
      r_idx  <- (N - n_times) + s
      w_d    <- compute_w_dhrp(sparse_cov_array[,,s], tau, s)
      turn   <- sum(abs(w_d - prev_d))
      
      if (turn > threshold_turn) {
        tx       <- calculate_transaction_costs(
          w_d, prev_d, prev_aum,
          price_mat[r_idx, ], vol_mat[r_idx, ],
          sqrt(diag(sparse_cov_array[,,s])),
          spread_vec
        )
        cost     <- tx$portfolio_cost
        val_aum  <- prev_aum - cost
        g        <- sum(w_d * returns_mat[r_idx, ])
        rets_val[i] <- g * (val_aum / prev_aum) - (cost / prev_aum)
        prev_d   <- w_d
        prev_aum <- val_aum * exp(g)
      } else {
        g           <- sum(prev_d * returns_mat[r_idx, ])
        rets_val[i] <- g
        prev_aum    <- prev_aum * exp(g)
        # prev_d remains unchanged
      }
    }
    
    ann_r <- mean(rets_val, na.rm = TRUE) * 252
    ann_v <- sd(rets_val, na.rm = TRUE)   * sqrt(252)
    sh    <- ann_r / ann_v
    
    if (!is.na(sh) && sh > best_sh) {
      best_sh  <- sh
      best_tau <- tau
    }
  }
  
  # --- Testing block using best_tau ---
  test_idx    <- t0:min(n_times, t0 + test_block - 1)
  prev_d_aum  <- aum_dhrp[t0 - 1]
  prev_n_aum  <- aum_naive[t0 - 1]
  prev_e_aum  <- aum_ew[t0 - 1]
  
  for (s in test_idx) {
    r_idx       <- (N - n_times) + s
    day_counter <- day_counter + 1
    
    # --- Equal‐weight (EW) strategy ---
    r_vec   <- returns_mat[r_idx, ]
    port_ret <- sum(w_ew_prev * r_vec)
    w_drift  <- (w_ew_prev * exp(r_vec)) / sum(w_ew_prev * exp(r_vec))
    w_target <- rep(1 / p_assets, p_assets)
    
    if (day_counter %% 25 == 0) {
      tx       <- calculate_transaction_costs(
        w_target, w_drift, prev_e_aum,
        price_mat[r_idx, ], vol_mat[r_idx, ],
        sqrt(diag(sparse_cov_array[,,s])),
        spread_vec
      )
      cost_ew[s] <- tx$portfolio_cost
      aac_ew     <- prev_e_aum - cost_ew[s]
      rets_ew[s] <- port_ret * (aac_ew / prev_e_aum) - (cost_ew[s] / prev_e_aum)
      aum_ew[s]  <- aac_ew * exp(port_ret)
      w_ew_prev  <- w_target
    } else {
      rets_ew[s] <- port_ret
      aum_ew[s]  <- prev_e_aum * exp(port_ret)
      w_ew_prev  <- w_drift
    }
    prev_e_aum <- aum_ew[s]
    
    # --- Naïve‐GMV strategy ---
    w_n   <- compute_w_naive(naive_sparse_cov_array[,,s])
    turnn <- sum(abs(w_n - w_naive_prev))
    
    if (turnn > threshold_turn) {
      tx          <- calculate_transaction_costs(
        w_n, w_naive_prev, prev_n_aum,
        price_mat[r_idx, ], vol_mat[r_idx, ],
        sqrt(diag(naive_sparse_cov_array[,,s])),
        spread_vec
      )
      cost_naive[s] <- tx$portfolio_cost
      aac           <- prev_n_aum - cost_naive[s]
      g             <- sum(w_n * returns_mat[r_idx, ])
      rets_naive[s] <- g * (aac / prev_n_aum) - (cost_naive[s] / prev_n_aum)
      aum_naive[s]  <- aac * exp(g)
      w_naive_prev  <- w_n
    } else {
      g              <- sum(w_naive_prev * returns_mat[r_idx, ])
      rets_naive[s] <- g
      aum_naive[s]  <- prev_n_aum * exp(g)
      # w_naive_prev unchanged
    }
    prev_n_aum <- aum_naive[s]
    
    # --- DHRP strategy (using sparse_cov_array and best_tau) ---
    w_d    <- compute_w_dhrp(sparse_cov_array[,,s], best_tau, s)
    turn_d <- sum(abs(w_d - w_dhrp_prev))
    
    if (turn_d > threshold_turn) {
      tx             <- calculate_transaction_costs(
        w_d, w_dhrp_prev, prev_d_aum,
        price_mat[r_idx, ], vol_mat[r_idx, ],
        sqrt(diag(sparse_cov_array[,,s])),
        spread_vec
      )
      cost_dhrp[s]   <- tx$portfolio_cost
      aac_dhrp       <- prev_d_aum - cost_dhrp[s]
      g              <- sum(w_d * returns_mat[r_idx, ])
      rets_dhrp[s]   <- g * (aac_dhrp / prev_d_aum) - (cost_dhrp[s] / prev_d_aum)
      aum_dhrp[s]    <- aac_dhrp * exp(g)
      w_dhrp_prev    <- w_d
    } else {
      g              <- sum(w_dhrp_prev * returns_mat[r_idx, ])
      rets_dhrp[s]   <- g
      aum_dhrp[s]    <- prev_d_aum * exp(g)
      # w_dhrp_prev updated to w_d (equivalent to no‐turn case)
      w_dhrp_prev    <- w_d
    }
    prev_d_aum <- aum_dhrp[s]
  }
}

# -----------------------------------
# 10. Results: compute performance stats
# -----------------------------------

idx         <- which(!is.na(rets_dhrp))
rets_xts    <- xts(
  cbind(
    DHRP  = rets_dhrp[idx],
    Naive = rets_naive[idx],
    EW    = rets_ew[idx]
  ),
  order.by = dates[(N - n_times) + idx]
)
simple_rets_xts <- exp(rets_xts) - 1

# Performance summary plot
charts.PerformanceSummary(simple_rets_xts,
                          main = "DHRP vs Naïve‐GMV vs EW (with Costs)")

# Tabulate metrics
df_perf <- data.frame(
  Strategy     = colnames(simple_rets_xts),
  Ann.Return   = round(as.numeric(Return.annualized( simple_rets_xts, 252)), 3),
  Ann.Vol      = round(as.numeric(StdDev.annualized(   simple_rets_xts, 252)), 3),
  Sharpe       = round(as.numeric(SharpeRatio.annualized(simple_rets_xts, 0, 252)), 3),
  Sortino      = round(as.numeric(SortinoRatio(            simple_rets_xts, MAR = 0)),    3),
  Max.Drawdown = round(as.numeric(maxDrawdown(             simple_rets_xts)),          3)
)
print(df_perf)

# Cumulative costs summary
cum_costs <- data.frame(
  Strategy            = c("DHRP", "Naive", "EW"),
  Total_Cost          = c(sum(cost_dhrp[idx], na.rm = TRUE),
                          sum(cost_naive[idx], na.rm = TRUE),
                          sum(cost_ew[idx],   na.rm = TRUE)),
  Avg_Cost_Per_Trade  = c(
    mean(cost_dhrp[idx][cost_dhrp[idx] > 0], na.rm = TRUE),
    mean(cost_naive[idx][cost_naive[idx] > 0], na.rm = TRUE),
    mean(cost_ew[idx][cost_ew[idx] > 0], na.rm = TRUE)
  ),
  Pct_of_Initial_AUM  = c(
    sum(cost_dhrp[idx], na.rm = TRUE) / initial_aum * 100,
    sum(cost_naive[idx], na.rm = TRUE) / initial_aum * 100,
    sum(cost_ew[idx], na.rm = TRUE)   / initial_aum * 100
  )
)
print(cum_costs)
