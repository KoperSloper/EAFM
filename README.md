# 📈 Sparse Covariance Portfolio Backtest with Transaction Costs (Including Market Impact)

This R script implements a backtesting framework for constructing and evaluating portfolios based on **sparse covariance matrix estimates** using **Realized GARCH** and **Graphical Lasso**. It compares the performance of three portfolio strategies—**Constrained Hierarchical Risk Parity (DHRP)**, **Naive Global Minimum Variance (Naive‐GMV)**, and **Equally Weighted (EW)**—while explicitly accounting for transaction costs (both bid‐ask spread and power‐law market impact).

---

## 🚀 Overview

1. **Data Processing**  
   - Filter out assets with missing/invalid returns or incomplete date coverage.  
   - Compute daily returns and realized volatility (Parkinson estimator).

2. **Covariance and Correlation Estimation**  
   - **RealGARCH (univariate)**: Model each asset’s realized volatility.  
   - **aDCC (multivariate)**: Estimate time‐varying correlation matrices using Dynamic Conditional Correlation.  
   - **Graphical Lasso (Glasso)**: Produce a sparse inverse covariance (precision) matrix from the regularized covariance estimate.

3. **Portfolio Construction**  
   - **DHRP (Constrained Hierarchical Risk Parity)**:  
     Build a hierarchical clustering tree based on the sparse precision matrix, then allocate risk in a “top‐down” fashion.  
   - **Naive‐GMV**:  
     Solve a Global Minimum Variance problem by inverting the (sparse) covariance matrix or precision matrix, under simple constraints (no shorting, unit sum).  
   - **EW (Equally Weighted)**:  
     Allocate one‐over‐p to each of p assets—benchmarked for comparison.

4. **Transaction Costs (Spread + Market Impact)**  
   - **Bid‐Ask Spread**:  
     Each share traded incurs half of the bid‐ask spread multiplied by the number of shares traded (i.e., volume).  
   - **Market Impact (Power‐Law)**:  
     When a trade is large relative to daily volume, it moves the price. We model this using a power‐law function of the dollar volume traded. A constant multiplier (α) and exponent (β > 1) govern how impact grows with trade size.  
   - **Combined Cost**:  
     For every portfolio rebalancing, we compute:  
     1. The spread cost, which equals half the bid‐ask spread times the trade volume.  
     2. The impact cost, which is α times the dollar volume raised to the power β.  
     The total transaction cost is the sum of spread cost plus impact cost. Whenever portfolio weights change (i.e., turnover exceeds a threshold), these costs are applied and subtracted from the portfolio’s assets under management (AUM).

5. **Backtesting / Performance Evaluation**  
   - Rolling estimation windows (e.g., 1,000 trading days), periodic rebalancing (e.g., every 25 days).  
   - Calculate **net returns** after deducting transaction costs each day.  
   - Compute and report:  
     - Annualized Return & Volatility  
     - Sharpe Ratio & Sortino Ratio  
     - Maximum Drawdown  
     - Cumulative Turnover & Total Costs  

---
