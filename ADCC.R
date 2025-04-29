library(tidyverse)
library(dbscan)
library(ggplot2)
library(zoo)
library(PerformanceAnalytics)
library(rmgarch)
library(covglasso)
library(knitr)
library(kableExtra)

data_stocks <- more_stocks

data_stocks_wide <- data_stocks %>%
  pivot_wider(
    names_from = PERMNO,     
    values_from = RET        
  )

data_stocks_wide <- data_stocks_wide %>%
  mutate(date = as.Date(date))

data_clean <- data_stocks_wide %>%
  select(
    date,
    where(~ !any(is.na(.)))
  )

data_clean <- data_clean %>%
  select(-`79996`)

data_logret_wide <- data_clean %>%
  mutate(across(-date, ~ as.numeric(.))) %>%
  mutate(across(-date, ~ log1p(.)))

returns_matrix <- as.matrix(data_logret_wide %>% select(-date))
rownames(returns_matrix) <- as.character(data_logret_wide$date)

univ_spec <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

multiv_spec <- dccspec(
  uspec = multispec(replicate(ncol(returns_matrix), univ_spec)),
  dccOrder = c(1, 1),
  model = "aDCC",
  distribution = "mvnorm"
)

# Test size is 1000 days and we use 1000 previous days for a prediction
rolling_fit <- dccroll(
  spec = multiv_spec,
  data = returns_matrix,
  n.ahead = 1,
  forecast.length = 1000,
  refit.every = 5,
  refit.window = "moving",
  window.size = 1000,
  fit.control = list(eval.se = FALSE)
)

sigma_mat <- sigma(rolling_fit)
cov_array <- rcov(rolling_fit)
corr_array <- rcor(rolling_fit)

# Making the covariance matrices sparser
p <- dim(cov_array)[1]
Tn <- dim(cov_array)[3]

sparse_cov_array  <- array(NA, dim = c(p, p, Tn))
sparse_corr_array <- array(NA, dim = c(p, p, Tn))

eps <- 1e-6
for (t in seq_len(Tn)) {
  S_t <- cov_array[,,t]
  # Added a tiny epsilon to the diagonal so S_t + eps·I is strictly PD
  S_reg <- S_t + diag(eps, ncol(S_t))  
  fit   <- covglasso(S = S_reg,
                     n      = 1000,
                     lambda = NULL,
                     path   = FALSE)
  sparse_cov_array[,,t]  <- fit$sigma
  sparse_corr_array[,,t] <- cov2cor(fit$sigma)
}

# See how many were "zeroed out" by the lasso penalty
count_zeros <- function(mat, exclude_diag = TRUE) {
  if (exclude_diag) {
    diag(mat) <- NA
  }
  sum(mat == 0, na.rm = TRUE)
}

zeros_raw   <- numeric(Tn)
zeros_sparse <- numeric(Tn)

for (t in seq_len(Tn)) {
  zeros_raw[t]    <- count_zeros(cov_array[,,t])
  zeros_sparse[t] <- count_zeros(sparse_cov_array[,,t])
}

mean_zeros_raw    <- mean(zeros_raw)
mean_zeros_sparse <- mean(zeros_sparse)

cat("Average number of zeros (raw):", mean_zeros_raw, "\n")
cat("Average number of zeros (sparse):", mean_zeros_sparse, "\n")

# So on average making the covariance matrices sparser resulted in 370 values being zeroed out

# parameters
lambda <- 0.5      
minPts <- 5        
eps_cl <- 0.003     

labels_mat <- matrix(NA_integer_, nrow = Tn, ncol = p)
colnames(labels_mat) <- colnames(returns_matrix)
rownames(labels_mat) <- paste0("t", seq_len(Tn))

for(t in seq_len(Tn)) {
  vol_t  <- as.numeric(sigma_mat[t, ])  
  corr_t <- sparse_corr_array[,, t]   
  
  vol_diff_mat  <- abs(outer(vol_t, vol_t, "-"))

  corr_abs_mat  <- abs(corr_t)
  
  D <- lambda * vol_diff_mat + (1 - lambda) * corr_abs_mat
  
  dist_obj <- as.dist(D)
  
  opt <- optics(dist_obj, minPts = minPts, eps = 0.05, search = "dist")
  
  labels_mat[t, ] <- extractDBSCAN(opt, eps_cl = eps_cl)$cluster
}

plot_vol_by_cluster <- function(t0,
                                sigma_mat,
                                labels_mat,
                                date_labels = rownames(labels_mat)) {
  df <- data.frame(
    Stock   = colnames(labels_mat),
    Vol     = as.numeric(sigma_mat[t0, ]),
    Cluster = factor(labels_mat[t0, ])
  )
  ggplot(df, aes(x = reorder(Stock, Vol), y = Vol, fill = Cluster)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(
      title    = sprintf("Volatilities & clusters (t = %d)", t0),
      subtitle = sprintf("Date: %s", date_labels[t0]),
      x        = NULL,
      y        = expression(sigma[i])
    )
}

plot_vol_by_cluster(t0 = 300, sigma_mat, labels_mat)

avg_within_cluster_corr2 <- function(t0,
                                     corr_array,
                                     labels_mat,
                                     include_noise = TRUE) {
  corr_t <- corr_array[,, t0]
  labels <- labels_mat[t0, ]
  
  if(include_noise) {
    clusters <- sort(unique(labels))    
  } else {
    clusters <- sort(unique(labels[labels > 0]))
  }
  
  result <- setNames(
    rep(NA_real_, length(clusters)),
    paste0("cluster_", clusters)
  )
  
  for(c in clusters) {
    members <- which(labels == c)
    if(length(members) >= 2) {
      subcorr <- corr_t[members, members]
      result[paste0("cluster_", c)] <- mean(subcorr[lower.tri(subcorr)])
    } 
  }
  
  result
}

avg_within_cluster_corr2(
  t0 = 297,
  corr_array = sparse_corr_array,
  labels_mat = labels_mat,
  include_noise = TRUE
)

# Build and compare portfolios
N             <- nrow(returns_matrix)
forecast_rows <- (N - Tn + 1):N
p             <- ncol(returns_matrix)

ret_strat <- numeric(Tn)
ret_naive <- numeric(Tn)

window.size <- 500

for (t in seq_len(Tn)) {
  actual_ret <- returns_matrix[ forecast_rows[t], ]
  
  # Strategy portfolio: lowest‐vol cluster, sparse covariance
  vols    <- sigma_mat[t, ]
  labs    <- labels_mat[t, ]
  avg_vol <- tapply(vols, labs, mean)
  best_cl <- as.integer(names(avg_vol)[which.min(avg_vol)])
  idx     <- which(labs == best_cl)
  
  S_sp   <- sparse_cov_array[idx, idx, t]
  invS   <- solve(S_sp)
  w_sp   <- invS %*% rep(1, length(idx))
  w_sp   <- as.numeric(w_sp / sum(w_sp))
  ret_strat[t] <- sum(w_sp * actual_ret[idx])
  
  #bNaïve GMV: all stocks, historical sample covariance
  hist_idx <- (forecast_rows[t] - window.size):(forecast_rows[t] - 1)
  S_hist   <- cov(returns_matrix[hist_idx, ])
  invH     <- solve(S_hist)
  w_naive  <- invH %*% rep(1, p)
  w_naive  <- as.numeric(w_naive / sum(w_naive))
  ret_naive[t] <- sum(w_naive * actual_ret)
}

N             <- nrow(returns_matrix)
Tn            <- dim(cov_array)[3]
forecast_rows <- (N - Tn + 1):N

forecast_dates <- as.Date(rownames(returns_matrix)[forecast_rows])

rets_xts <- xts(
  x        = cbind(Strategy = ret_strat, Naive = ret_naive),
  order.by = forecast_dates
)

charts.PerformanceSummary(
  rets_xts,
  main = "Lowest-Vol Cluster Strategy vs. Naïve GMV"
)

annRet  <- Return.annualized(rets_xts, scale = 252)           
annVol  <- StdDev.annualized(rets_xts,  scale = 252)           

sharpe  <- SharpeRatio.annualized(rets_xts, Rf = 0, scale = 252)

sortino <- SortinoRatio(rets_xts, MAR = 0)

mdd     <- maxDrawdown(rets_xts)

perf_df <- data.frame(
  Strategy     = c("Cluster-GMV", "Naïve-GMV"),
  Ann.Return   = round(as.numeric(annRet),    3),
  Ann.Vol      = round(as.numeric(annVol),    3),
  Sharpe       = round(as.numeric(sharpe),    3),
  Sortino      = round(as.numeric(sortino),   3),
  Max.Drawdown = round(as.numeric(mdd),       3)
)

print(perf_df)
