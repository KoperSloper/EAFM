library(tidyverse)
library(xts)
library(rugarch)
library(ggplot2)
library(zoo)
library(PerformanceAnalytics)

data_stocks <- new_stocks_daily

# Make it wide format
data_stocks_wide <- data_stocks %>%
  pivot_wider(
    names_from = PERMNO,     
    values_from = RET        
  )

# Correct formatting for date
data_stocks_wide <- data_stocks_wide %>%
  mutate(date = as.Date(date))

# Remove stocks where we have any missing values
data_clean <- data_stocks_wide %>%
  select(
    date,
    where(~ !any(is.na(.)))
  )

# Make it log returns
data_logret_wide <- data_clean %>%
  mutate(across(-date, ~ log1p(.)))

# Define EGARCH spec once outside the function
spec_egarch <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(1,1)),
  mean.model         = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "norm"
)

# Create a function to process each stock
process_stock <- function(stock_name, data_wide, spec) {
  # Convert to xts
  ts_data <- xts(
    data_wide[[stock_name]], 
    order.by = as.Date(data_wide$date)
  )
  
  # Run ugarchroll
  roll_result <- ugarchroll(
    spec            = spec,
    data            = ts_data,
    n.ahead         = 1,
    forecast.length = 250,
    refit.every     = 1,
    refit.window    = "recursive"
  )
  
  return(roll_result)
}

stock_names <- setdiff(names(data_logret_wide), "date")

results_list <- lapply(stock_names, process_stock, data_wide = data_logret_wide, spec = spec_egarch)

names(results_list) <- stock_names

volatility_df <- data.frame()

for (stock_name in names(results_list)) {
  forecast_data <- as.data.frame(results_list[[stock_name]], which="density")
  
  stock_df <- data.frame(
    date = as.Date(rownames(forecast_data)),
    stock = stock_name,
    sigma = forecast_data$Sigma
  )
  
  volatility_df <- rbind(volatility_df, stock_df)
}

volatility_df <- volatility_df %>% 
  dplyr::arrange(volatility_df$date, volatility_df$stock)

volatility_wide <- volatility_df %>%
  pivot_wider(
    names_from = stock,
    values_from = sigma,
    id_cols = date
  )


apply_clustering_by_date <- function(volatility_wide, 
                                     minPts = 5,  
                                     eps    = 10,    
                                     eps_cl = 0.5) {
  
  volatility_long <- volatility_wide %>%
    pivot_longer(
      cols = -date,
      names_to  = "stock",
      values_to = "sigma"
    ) %>%
    mutate(date = as.Date(date))

  dates <- sort(unique(volatility_long$date))

  cluster_results  <- tibble(
    date            = as.Date(character()),
    min_vol_cluster = integer(),
    num_stocks      = integer(),
    avg_volatility  = numeric()
  )
  portfolio_stocks <- list()

  for (d in dates) {
    vol_data <- volatility_long %>%
      filter(date == d, !is.na(sigma))
    
    if (nrow(vol_data) >= minPts) {
      optics_res <- optics(
        matrix(vol_data$sigma, ncol = 1),
        eps    = eps,
        minPts = minPts
      )
      clusters <- extractDBSCAN(optics_res, eps_cl = eps_cl)$cluster
      vol_data$cluster <- clusters

      cluster_summary <- vol_data %>%
        group_by(cluster) %>%
        summarize(
          avg_sigma = mean(sigma),
          count     = n(),
          .groups   = "drop"
        ) %>%
        arrange(avg_sigma)

      valid <- cluster_summary %>% filter(cluster > 0)
      if (nrow(valid) > 0) {
        chosen_cluster <- valid$cluster[1]
        sel_stocks <- vol_data %>%
          filter(cluster == chosen_cluster) %>%
          pull(stock)
      } else {
        sel_stocks <- vol_data %>%
          filter(cluster != 0) %>%
          pull(stock)
        if (length(sel_stocks) == 0) {
          sel_stocks <- vol_data %>%
            arrange(sigma) %>%
            slice_head(prop = 0.1) %>%
            pull(stock)
        }
        chosen_cluster <- NA
      }

      portfolio_stocks[[ as.character(d) ]] <- sel_stocks

      cluster_results <- bind_rows(
        cluster_results,
        tibble(
          date            = as.Date(d),    # ← change: ensure this is Date
          min_vol_cluster = ifelse(is.na(chosen_cluster), -1L, chosen_cluster),
          num_stocks      = length(sel_stocks),
          avg_volatility  = mean(vol_data$sigma[vol_data$stock %in% sel_stocks])
        )
      )
    }
  }
  
  names(portfolio_stocks) <- as.character(cluster_results$date)
  
  list(
    cluster_results  = cluster_results,
    portfolio_stocks = portfolio_stocks
  )
}

clustering_results <- apply_clustering_by_date(
  volatility_wide,
  minPts = 2,
  eps    = 0.005,
  eps_cl = 0.04
)

visualize_volatility_selection <- function(volatility_wide, clustering_results, date_to_viz) {
  volatility_long <- volatility_wide %>%
    pivot_longer(cols = -date, 
                 names_to = "stock", 
                 values_to = "sigma")
  
  vol_data <- volatility_long %>% 
    filter(date == date_to_viz) %>%
    filter(!is.na(sigma))
  
  if (!as.character(date_to_viz) %in% names(clustering_results$portfolio_stocks)) {
    stop(paste("No clustering results available for date:", date_to_viz))
  }
  
  selected_stocks <- clustering_results$portfolio_stocks[[as.character(date_to_viz)]]
  
  vol_data <- vol_data %>%
    mutate(in_portfolio = stock %in% selected_stocks)

  vol_data$group <- ifelse(vol_data$in_portfolio, "Selected Portfolio", "Other Stocks")

  p <- ggplot(vol_data, aes(x = reorder(stock, sigma), y = sigma, color = group)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Other Stocks" = "#FFB6C1", "Selected Portfolio" = "steelblue")) +
    labs(
      title = paste("Stock Volatility Distribution on", date_to_viz),
      subtitle = paste(sum(vol_data$in_portfolio), "stocks in low-volatility portfolio out of", nrow(vol_data), "total"),
      x = "Stocks (ordered by volatility)",
      y = "Volatility (sigma)",
      color = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "bottom"
    )
  
  print(p)

  return(invisible(vol_data))
}

dot_plot_data <- visualize_volatility_selection(
  volatility_wide, 
  clustering_results, 
  "2024-08-06"
)

compare_portfolio_strategies <- function(data_logret_wide, clustering_results) {
  all_dates <- sort(unique(data_logret_wide$date))
  
  overlap_dates <- as.character(all_dates)[as.character(all_dates) %in% names(clustering_results$portfolio_stocks)]
  cat("Number of overlapping dates:", length(overlap_dates), "\n")
  if (length(overlap_dates) > 0) {
    cat("First few overlapping dates:", head(overlap_dates), "\n")
  } else {
    cat("NO OVERLAPPING DATES FOUND!\n")
  }

  cluster_portfolio_returns <- data.frame(date = all_dates[-1], return = NA_real_)
  all_stocks_returns       <- data.frame(date = all_dates[-1], return = NA_real_)

  portfolio_compositions <- data.frame()
  successful_counts      <- 0
  
  for (i in seq_along(all_dates)[-length(all_dates)]) {
    current_date <- all_dates[i]
    next_date    <- all_dates[i + 1]

    date_key <- as.character(current_date)
    
    if (! date_key %in% names(clustering_results$portfolio_stocks)) next
    
    selected_stocks <- clustering_results$portfolio_stocks[[date_key]]
    selected_stocks <- as.character(selected_stocks)  # if factor

    curr_ret <- data_logret_wide %>% filter(date == current_date)
    next_ret <- data_logret_wide %>% filter(date == next_date)
    
    if (nrow(curr_ret) && nrow(next_ret)) {
      valid_sel <- intersect(selected_stocks, names(curr_ret)[-1])
      if (length(valid_sel)) {
        portfolio_compositions <- rbind(
          portfolio_compositions,
          data.frame(
            date                = current_date,
            next_date           = next_date,
            num_selected_stocks = length(valid_sel),
            num_all_stocks      = ncol(curr_ret) - 1
          )
        )
        
        cluster_portfolio_returns$return[cluster_portfolio_returns$date == next_date] <-
          next_ret %>% select(all_of(valid_sel)) %>% unlist() %>% mean(na.rm = TRUE)
        
        all_stocks_returns$return[all_stocks_returns$date == next_date] <-
          next_ret %>% select(-date) %>% unlist() %>% mean(na.rm = TRUE)
        
        successful_counts <- successful_counts + 1
      }
    }
  }
  
  cat("Successfully formed portfolios for", successful_counts, "dates\n")

  cluster_portfolio_returns <- cluster_portfolio_returns %>% filter(!is.na(return))
  all_stocks_returns        <- all_stocks_returns       %>% filter(!is.na(return))

  combined_returns <- inner_join(
    cluster_portfolio_returns  %>% rename(cluster_return    = return),
    all_stocks_returns         %>% rename(all_stocks_return = return),
    by = "date"
  )
  cat("Number of dates with combined returns:", nrow(combined_returns), "\n")

  combined_returns <- combined_returns %>%
    mutate(
      cum_cluster_return     = cumprod(1 + cluster_return)     - 1,
      cum_all_stocks_return  = cumprod(1 + all_stocks_return)  - 1
    )

  portfolio_metrics <- tibble(
    strategy             = c("Clustering Low-Vol", "All Stocks EW"),
    annualized_return    = c(
      mean(combined_returns$cluster_return,    na.rm = TRUE) * 252,
      mean(combined_returns$all_stocks_return, na.rm = TRUE) * 252
    ),
    annualized_volatility = c(
      sd(combined_returns$cluster_return,    na.rm = TRUE) * sqrt(252),
      sd(combined_returns$all_stocks_return, na.rm = TRUE) * sqrt(252)
    )
  ) %>%
    mutate(sharpe_ratio = annualized_return / annualized_volatility)

  plot_data <- combined_returns %>%
    select(date, cum_cluster_return, cum_all_stocks_return) %>%
    pivot_longer(
      starts_with("cum_"),
      names_to  = "strategy",
      values_to = "cumulative_return"
    ) %>%
    mutate(strategy = recode(
      strategy,
      cum_cluster_return    = "Clustering Low-Vol",
      cum_all_stocks_return = "All Stocks EW"
    ))
  
  p <- ggplot(plot_data, aes(date, cumulative_return, color = strategy)) +
    geom_line(linewidth = 1) +
    labs(
      title = "Comparison of Portfolio Strategies: Cumulative Returns",
      x     = "Date",
      y     = "Cumulative Return",
      color = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  print(portfolio_metrics)
  
  list(
    daily_returns          = combined_returns,
    portfolio_metrics      = portfolio_metrics,
    portfolio_compositions = portfolio_compositions,
    plot                   = p
  )
}


comparison_results <- compare_portfolio_strategies(
  data_logret_wide,
  clustering_results
)

