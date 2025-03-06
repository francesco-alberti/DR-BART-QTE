library(ggplot2)   
library(splines)  
library(quantreg)  
library(drbart)

###################### ILLUSTRATION 1: 2 programs ATE vs QTE ###################
set.seed(42)
# Weibull parameters
lambda_0 <- 2  
k_0 <- 1.5
# Closed-form CDFs
cdf_weibull <- function(x, lambda, k) {
  1 - exp(- (x / lambda)^k)
}
# Closed-form inverse CDFs (quantile functions)
quantile_weibull <- function(p, lambda, k) {
  lambda * (-log(1 - p))^(1 / k)
}
# Quantile sequence
tau_seq <- seq(0.01, 0.99, length.out = 100)
# Generate Y_0 quantiles
Q_Y0 <- quantile_weibull(tau_seq, lambda_0, k_0)
# Program A: parameters and Y_1 quantiles
lambda_A <- lambda_0 * 1.2  
k_A <- k_0 * 1.3  
Q_Y1_A <- quantile_weibull(tau_seq, lambda_A, k_A)
# Program B: parameters and Y_1 quantiles
lambda_B <- lambda_0 * 1.2  
k_B <- k_0 * 0.7 
Q_Y1_B <- quantile_weibull(tau_seq, lambda_B, k_B)
# True QTEs
QTE_A_true <- Q_Y1_A - Q_Y0
QTE_B_true <- Q_Y1_B - Q_Y0
# True ATEs
ATE_A_true <- mean(QTE_A_true)
ATE_B_true <- mean(QTE_B_true)
# Sequence of x values for plotting true CDFs
x_vals <- seq(0, 10, length.out = 500)
# True CDFs
cdf_Y0 <- cdf_weibull(x_vals, lambda_0, k_0)
cdf_Y1_A <- cdf_weibull(x_vals, lambda_A, k_A)
cdf_Y1_B <- cdf_weibull(x_vals, lambda_B, k_B)
# Plot true CDFs for Program A and Program B, separate plots
df_CDF_A <- data.frame(X = x_vals, Y0 = cdf_Y0, Y1_A = cdf_Y1_A)
df_CDF_B <- data.frame(X = x_vals, Y0 = cdf_Y0, Y1_B = cdf_Y1_B)

plot_CDF_A <- ggplot(df_CDF_A, aes(x = X)) +
  geom_line(aes(y = Y0, color = "Control"), size = 1) +
  geom_line(aes(y = Y1_A, color = "Treatment"), size = 1) +
  labs(title = "CDFs for Program A",
       x = "Wage (Y)", y = "CDF", color = "Group") +
  theme_minimal()

plot_CDF_B <- ggplot(df_CDF_B, aes(x = X)) +
  geom_line(aes(y = Y0, color = "Control"), size = 1) +
  geom_line(aes(y = Y1_B, color = "Treatment"), size = 1) +
  labs(title = "CDFs for Program B",
       x = "Wage (Y)", y = "CDF", color = "Group") +
  theme_minimal()

# Plot True QTE and ATE for both programs A and B, same plot
df_QTE <- data.frame(
  Quantile = tau_seq,
  QTE_A = QTE_A_true,
  QTE_B = QTE_B_true
)

plot_QTE_ATE <- ggplot(df_QTE, aes(x = Quantile)) +
  geom_line(aes(y = QTE_A, color = "A"), size = 1) +  
  geom_line(aes(y = QTE_B, color = "B"), size = 1) +  
  geom_hline(yintercept = ATE_A_true, linetype = "dashed", color = "blue", size = 0.5) + 
  geom_hline(yintercept = ATE_B_true, linetype = "dashed", color = "red", size = 0.5) +
  labs(title = "", #QTEs and ATEs for Programs A and B
       x = "Quantile", y = "Treatment Effect",
       color = "Program") +  
  scale_color_manual(values = c("blue", "red"))+
  theme_minimal()+ 
  theme(
    legend.position = c(0.1, 0.8),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  ) 
print(plot_CDF_A)  
print(plot_CDF_B)  
print(plot_QTE_ATE) 

###################### ILLUSTRATION 2: quantile crossing - quantreg vs dr-bart ###################
set.seed(999)
# Dataset
n <- 50
x <- sort(runif(n, 0, 10))  
y <- 3 + sin(2*x) - 0.2*x^2 + rnorm(n, sd = 5)

#------ QUANTILE REGRESSION ------
# Quantiles of interest
taus <- c(0.4, 0.5, 0.6)
# Quantile regression model with B-Splines and df=5
models <- lapply(taus, function(tau) rq(y ~ bs(x, df = 5), tau = tau))  
x_pred <- seq(min(x), max(x), length.out = 100) 
preds_qr <- sapply(models, function(model) predict(model, newdata = data.frame(x = x_pred)))
# Plot conditional quantiles estimated via B-Spline quantile regression
pred_df_qr <- data.frame(x = rep(x_pred, times = length(taus)),
                         y = as.vector(preds_qr),
                         quantile = rep(paste0(taus * 100, "%"), each = length(x_pred)))

plot_quantreg <- ggplot() +
  geom_point(aes(x, y), alpha = 0.3) +   
  geom_line(data = pred_df_qr, aes(x, y, color = quantile), linewidth = 0.7) +  
  labs(title = "", # Conditional Quantiles via B-Splines Quantile Regression
       x = "X",
       y = "Y",
       color = "Quantile") +
  theme_minimal() +
  ylim(min(y), max(y)) + 
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   
print(plot_quantreg)
                   
#------ DR-BART ------
fit <- drbart(y, 
              x,
              nburn = 2000,
              nsim = 2000,
              variance = "const")

xpred <- matrix(seq(min(x), max(x), length.out = 20))
ygrid <- seq(min(y), max(y), length.out = 100)

preds <- predict(fit, 
                 xpred, 
                 ygrid,
                 type = "distribution")
# Function to invert the CDF to obtain quantiles
get_quantile <- function(cdf, grid, probs) {
  sapply(probs, function(p) grid[which.min(abs(cdf - p))])
}
# Compute posterior mean of CDF estimates
cdf_mean <- apply(preds$preds, c(1,2), mean)  
# Compute conditional quantiles 
quantiles<- apply(cdf_mean, 1, function(cdf) get_quantile(cdf, ygrid, taus))
quantiles <- t(quantiles) 
# Apply smoothing spline to each quantile curve
smoothed_quantiles <- apply(quantiles, 2, function(q) {
  smooth.spline(xpred, q, spar = 0.35)$y  # Adjust `spar` for smoothness
})
# Plot smoothed conditional quantiles estimated via dr-bart
smoothed_pred_df <- data.frame(
  x = rep(xpred, times = length(taus)),
  y = as.vector(smoothed_quantiles),
  quantile = rep(paste0(taus * 100, "%"), each = length(xpred))
)

plot_smooth_drbart <- ggplot() +
  geom_point(aes(x, y), alpha = 0.3) +  
  geom_line(data = smoothed_pred_df, aes(x = x, y = y, color = quantile), linewidth = 0.7) +  
  labs(title = "", #Smoothed Conditional Quantiles via DR-BART
       x = "X",
       y = "Y",
       color = "Quantile") +
  theme_minimal() +
  ylim(min(y), max(y)) + 
  theme(
    legend.position = c(0.8, 0.8),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5) 
  ) 
print(plot_smooth_drbart)

