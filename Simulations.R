library(drbart)     
library(ggplot2)     
library(qte)         
library(tidyverse)  
library(dplyr)       
library(reshape2)         
library(weights)

###################### SIMULATION 1: randomized experiment DR-BART QTE Gamma ###################
#------ TRUE CDFs and QTE ------ 
# Gamma distribution parameters
shape_0 <- 2  
scale_0 <- 1  
shape_1 <- 2  
scale_1 <- 1.5
# Sequence of Y values for CDFs
y_vals <- seq(0, 14, length.out = 500)

# Plot True CDFs
CDFy0_true <- pgamma(y_vals, shape = shape_0, scale = scale_0)  
CDFy1_true <- pgamma(y_vals, shape = shape_1, scale = scale_1) 

df_CDF_true <- data.frame(X = y_vals, Y0 = CDFy0_true, Y1 = CDFy1_true)

plot_CDF_true1 <- ggplot(df_CDF_true, aes(x = X)) +
  geom_line(aes(y = Y0, color = "Control"), size = 1) +
  geom_line(aes(y = Y1, color = "Treatment"), size = 1) +
  labs(title = "", #True CDF of Y
       x = "Y", y = "Probability", color = "Group") +
  theme_minimal()+ 
  theme(
    legend.position = c(0.80, 0.2),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  ) 

# Quantile sequence
tau_seq <- seq(0.01, 0.99, by = 0.01)
# True QTE using the inverse CDF (quantile function)
Q_Y0_true <- qgamma(tau_seq, shape = shape_0, scale = scale_0)  
Q_Y1_true <- qgamma(tau_seq, shape = shape_1, scale = scale_1)  
QTE_true <- Q_Y1_true - Q_Y0_true

# Plot True QTE
df_QTE_true <- data.frame(Quantile = tau_seq, QTE = QTE_true)

plot_QTE_true1 <- ggplot(df_QTE_true, aes(x = Quantile, y = QTE)) +
  geom_line(color = "blue", size = 1) +  
  labs(title = "", # True QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal() 

print(plot_CDF_true1)  
print(plot_QTE_true1) 

#------ EMPIRICAL CDFs and QTE ------ 
set.seed(42)
n <- 1000
# Treatment assignment  D (binary, randomized)
D <- rbinom(n, 1, 0.5)
# Potential outcomes, Y_0 and Y_1, independent of treatment D
Y_0 <- rgamma(n, shape = shape_0, scale = scale_0)   
Y_1 <- rgamma(n, shape = shape_1, scale = scale_1) 
# Observed outcome Y based on treatment assignment
Y <- ifelse(D == 1, Y_1, Y_0)
# Data frame
df <- data.frame(D = D, Y = Y)

# Plot empirical CDFs 
custom_colors <- c("Control" = "#F8766D", "Treatment" = "#00BFC4")
plot_CDF_empirical1 <- ggplot(df, aes(x = Y, color = factor(D, levels = c(0,1), labels = c("Control", "Treatment")))) +
  stat_ecdf(geom = "step", size = 1) +
  labs(title = "", x = "Y", y = "Probability", color = "Group") + #Empirical CDF of Y  
  scale_color_manual(values = custom_colors) +  
  theme_minimal()+ 
  theme(
    legend.position = c(0.80, 0.2),  # Legend in bottom-right corner inside plot
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  ) 

# Empirical quantiles for treated and control groups
Q_Y1_empirical <- quantile(Y[D == 1], probs = tau_seq, na.rm = TRUE)  
Q_Y0_empirical <- quantile(Y[D == 0], probs = tau_seq, na.rm = TRUE)  
# Empirical QTE 
QTE_empirical <- Q_Y1_empirical - Q_Y0_empirical

# Plot Empirical QTE
df_QTE_empirical <- data.frame(tau = tau_seq, QTE = QTE_empirical)

plot_QTE_empirical1 <- ggplot(df_QTE_empirical, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +                 
  labs(title = "",   #Empirical QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_CDF_empirical1)
print(plot_QTE_empirical1)

#------ DR-BART CDFs and QTE with CI ------ 
fit_DRBART <- drbart(Y,
                     D,
                     nburn = 2000,
                     nsim = 2000,
                     variance = "ux"
)
grid = seq(from = min(Y), to = max(Y), by = 0.05)
x_pred <- rbind(1,0)
preds_cdf <- predict(fit_DRBART, xpred = x_pred, ygrid = grid, type = "distribution")

# Extract mean CDF and 95% cred int 
CDFy_drbart <- apply(preds_cdf$preds, c(1,2), mean)
CDFy_drbart_2.5 <- apply(preds_cdf$preds, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
CDFy_drbart_97.5 <- apply(preds_cdf$preds, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
# Mean CDF for Treatment and Control
CDFy1_drbart <- CDFy_drbart[1, ]
CDFy0_drbart <- CDFy_drbart[2, ]
# 2.5  CDF for Treatment and Control
CDFy1_drbart_2.5 <- CDFy_drbart_2.5[1, ]
CDFy0_drbart_2.5 <- CDFy_drbart_2.5[2, ]
# 97.5  CDF for Treatment and Control
CDFy1_drbart_97.5 <- CDFy_drbart_97.5[1, ]
CDFy0_drbart_97.5 <- CDFy_drbart_97.5[2, ]

# Plot DR-BART CDFs with cred int
group_colors <- c("Control" = "#F8766D", "Treatment" = "#00BFC4")
df_CDF_drbart <- data.frame(
  grid,                 
  cdf = c(CDFy1_drbart, CDFy0_drbart),  
  cdf_2.5 = c(CDFy1_drbart_2.5, CDFy0_drbart_2.5),
  cdf_97.5 = c(CDFy1_drbart_97.5, CDFy0_drbart_97.5),
  group = rep(c("Treatment", "Control"), each = length(grid))
)

plot_CDF_drbart1 <- ggplot(df_CDF_drbart, aes(x = grid, y = cdf, color = group, fill = group)) +
  geom_ribbon(aes(ymin = cdf_2.5, ymax = cdf_97.5), alpha = 0.2, color = NA) +  
  geom_line(size = 1) +  
  scale_color_manual(values = group_colors) +  
  scale_fill_manual(values = group_colors) +  
  labs(
    title = "", #Estimated CDFs via DR-BART
    x = "Y",
    y = "Probability",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() + 
  theme(
    legend.position = c(0.8, 0.2),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )

print(plot_CDF_drbart1)
 
# Separate posterior samples for Y_1 and Y_0
preds1 <- preds_cdf$preds[1, , ]  
preds0 <- preds_cdf$preds[2, , ]
# Function to invert the CDF to obtain quantiles
get_quantile <- function(cdf, grid, probs) {
  sapply(probs, function(p) grid[which.min(abs(cdf - p))])  # Find closest CDF value and return corresponding Y
}
# Quantiles for each posterior draw
quant0 <- apply(preds0, 2, get_quantile, grid = grid, probs = tau_seq)  
quant1 <- apply(preds1, 2, get_quantile, grid = grid, probs = tau_seq)
# QTE for each posterior sample
qte_samples <- quant1 - quant0 
# Mean QTE
QTE_drbart <- apply(qte_samples, 1, mean)
# 95% cred int (2.5th and 97.5th percentiles)
QTE_drbart_2.5 <- apply(qte_samples, 1, quantile, probs = 0.025)
QTE_drbart_97.5 <- apply(qte_samples, 1, quantile, probs = 0.975)

# Plot DR-BART QTE with cred int
qte_drbart_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_drbart,
  QTE_2.5 = QTE_drbart_2.5,  
  QTE_97.5 = QTE_drbart_97.5 
)

plot_QTE_drbart1 <- ggplot(qte_drbart_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  
  labs(title = "", #DR-BART QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_drbart1)

#------ Firpo QTE with CI ------- 
Firpo <- ci.qte(Y ~ D,
                data=df,
                probs=tau_seq,
                se= TRUE,
                iters=100)
QTE_Firpo <- Firpo$qte
QTE_Firpo_2.5 <- Firpo$qte.lower
QTE_Firpo_97.5 <- Firpo$qte.upper 

# Plot Firpo QTE with conf int
qte_Firpo_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_Firpo,
  QTE_2.5 = QTE_Firpo_2.5,  
  QTE_97.5 = QTE_Firpo_97.5 
)

plot_QTE_Firpo1 <- ggplot(qte_Firpo_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  
  labs(title = "", #Firpo (2007) QTE 
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_Firpo1)

#------ TRUE QTE and Firpo QTE with CI ------
# Plot True QTE and Firpo QTE with conf int
plot_combined1 <- ggplot() +
  geom_ribbon(data = qte_Firpo_df, aes(x = tau, ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2, show.legend = FALSE) +
  geom_line(data = qte_Firpo_df, aes(x = tau, y = QTE, color = "Firpo (2007)"), size = 1) +
  geom_line(data = df_QTE_true, aes(x = Quantile, y = QTE, color = "True"), size = 0.7, linetype = "longdash") +
  labs(title = "", x = "Quantile", y = "Treatment Effect", color = "QTE") +
  scale_color_manual(values = c("True" = "red", "Firpo (2007)" = "blue")) +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_combined1)

#------ TRUE QTE and DR-BART QTE with CI ------
# Plot True QTE and DR-BART QTE with cred int
plot_combined2 <- ggplot() +
  geom_ribbon(data = qte_drbart_df, aes(x = tau, ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2, show.legend = FALSE) +
  geom_line(data = qte_drbart_df, aes(x = tau, y = QTE, color = "DR-BART"), size = 1) +
  geom_line(data = df_QTE_true, aes(x = Quantile, y = QTE, color = "True"), size = 0.7, linetype = "longdash") +
  labs(title = "", x = "Quantile", y = "Treatment Effect", color = "QTE") +
  scale_color_manual(values = c("True" = "red", "DR-BART" = "blue")) +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_combined2)

#------ RESULTS ------
# Distance measures for Empirical estimation
L1_distance_empirical <- sum(abs(QTE_true - QTE_empirical))
L2_distance_empirical <- sqrt(sum((QTE_true - QTE_empirical)^2))
Wasserstein_distance_empirical <- sum(abs(QTE_true - QTE_empirical)) / length(QTE_true)
max_abs_deviation_empirical <- max(abs(QTE_true - QTE_empirical))
ISE_empirical <- sum((QTE_true - QTE_empirical)^2) / length(QTE_true)

# Distance measures for DR-BART estimation
L1_distance_drbart <- sum(abs(QTE_true - QTE_drbart))
L2_distance_drbart <- sqrt(sum((QTE_true - QTE_drbart)^2))
Wasserstein_distance_drbart <- sum(abs(QTE_true - QTE_drbart)) / length(QTE_true)
max_abs_deviation_drbart <- max(abs(QTE_true - QTE_drbart))
ISE_drbart <- sum((QTE_true - QTE_drbart)^2) / length(QTE_true)

# Distance measures for Firpo estimation
L1_distance_Firpo <- sum(abs(QTE_true - QTE_Firpo))
L2_distance_Firpo <- sqrt(sum((QTE_true - QTE_Firpo)^2))
Wasserstein_distance_Firpo <- sum(abs(QTE_true - QTE_Firpo)) / length(QTE_true)
max_abs_deviation_Firpo <- max(abs(QTE_true - QTE_Firpo))
ISE_Firpo <- sum((QTE_true - QTE_Firpo)^2) / length(QTE_true)

# Results  
cat("L1 Distance DRBART:", L1_distance_drbart, "\n")
cat("L1 Distance EMPIRICAL:", L1_distance_empirical, "\n")
cat("L1 Distance Firpo:", L1_distance_Firpo, "\n")

cat("L2 Distance DRBART:", L2_distance_drbart, "\n")
cat("L2 Distance EMPIRICAL:", L2_distance_empirical, "\n")
cat("L2 Distance Firpo:", L2_distance_Firpo, "\n")

cat("Wasserstein Distance DRBART:", Wasserstein_distance_drbart, "\n")
cat("Wasserstein Distance EMPIRICAL:", Wasserstein_distance_empirical, "\n")
cat("Wasserstein Distance Firpo:", Wasserstein_distance_Firpo, "\n")

cat("Max Absolute Deviation DRBART:", max_abs_deviation_drbart, "\n")
cat("Max Absolute Deviation EMPIRICAL:", max_abs_deviation_empirical, "\n")
cat("Max Absolute Deviation Firpo:", max_abs_deviation_Firpo, "\n")

cat("ISE DRBART:", ISE_drbart, "\n\n")
cat("ISE EMPIRICAL:", ISE_empirical, "\n")
cat("ISE Firpo:", ISE_Firpo, "\n")
###################### SIMULATION 2: randomized experiment DR-BART QTE Gamma with additive and multiplicative noise ###################
#------ TRUE CDFs and QTE ------ 
# Gamma distribution parameters
shape_0 <- 2  
scale_0 <- 1  
shape_1 <- 2  
scale_1 <- 1.5
# Noise parameters
sigma_additive <- 2 
delta_multiplicative <- 1.5  
# Sequence of Y values for CDFs
y_vals <- seq(-5, 25, length.out = 500)
# Monte Carlo Approximation of True CDFs with noise
MC_samples <- 100000
Gamma_Y0 <- rgamma(MC_samples, shape = shape_0, scale = scale_0)
Gamma_Y1 <- rgamma(MC_samples, shape = shape_1, scale = scale_1)
Y0_samples <- Gamma_Y0 * runif(MC_samples, 1 - delta_multiplicative, 1 + delta_multiplicative) + rnorm(MC_samples, 0, sigma_additive)
Y1_samples <- Gamma_Y1 * runif(MC_samples, 1 - delta_multiplicative, 1 + delta_multiplicative) + rnorm(MC_samples, 0, sigma_additive)
CDFy0_true <- ecdf(Y0_samples)(y_vals)
CDFy1_true <- ecdf(Y1_samples)(y_vals)

# Plot True CDFs
df_CDF_true <- data.frame(X = y_vals, Y0 = CDFy0_true, Y1 = CDFy1_true)

plot_CDF_true2 <- ggplot(df_CDF_true, aes(x = X)) +
  geom_line(aes(y = Y0, color = "Control"), size = 1) +
  geom_line(aes(y = Y1, color = "Treatment"), size = 1) +
  labs(title = "", #True CDF of Y
       x = "Y", y = "Probability", color = "Group") +
  theme_minimal()

# Quantile sequence
tau_seq <- seq(0.01, 0.99, by = 0.01)
# True QTE using quantile function
Q_Y0_true <- quantile(Y0_samples, tau_seq)  
Q_Y1_true <- quantile(Y1_samples, tau_seq)  
QTE_true <- Q_Y1_true - Q_Y0_true

# Plot True QTE
df_QTE_true <- data.frame(Quantile = tau_seq, QTE = QTE_true)

plot_QTE_true2 <- ggplot(df_QTE_true, aes(x = Quantile, y = QTE)) +
  geom_line(color = "blue", size = 1) +  
  labs(title = "",  #True QTE and ATE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_CDF_true2)  
print(plot_QTE_true2)

#------ EMPIRICAL CDFs and QTE ------ 
set.seed(45)
n <- 1000
# Treatment assignment (randomized)
D <- rbinom(n, 1, 0.5)
# Potential outcomes with both additive & multiplicative noise
Y_0 <- rgamma(n, shape = shape_0, scale = scale_0) * runif(n, 1 - delta_multiplicative, 1 + delta_multiplicative) + rnorm(n, 0, sigma_additive)
Y_1 <- rgamma(n, shape = shape_1, scale = scale_1) * runif(n, 1 - delta_multiplicative, 1 + delta_multiplicative) + rnorm(n, 0, sigma_additive)
# Observed outcome based on treatment assignment
Y <- ifelse(D == 1, Y_1, Y_0)
df <- data.frame(D = D, Y = Y)

# Plot empirical CDFs 
custom_colors <- c("Control" = "#F8766D", "Treatment" = "#00BFC4")
plot_CDF_empirical2 <- ggplot(df, aes(x = Y, color = factor(D, levels = c(0,1), labels = c("Control", "Treatment")))) +
  stat_ecdf(geom = "step", size = 1) +
  labs(title = "", x = "Y", y = "Probability", color = "Group") +  #Empirical CDF of Y 
  scale_color_manual(values = custom_colors) +  
  theme_minimal()

# Empirical quantiles for treated and control groups
Q_Y1_empirical <- quantile(Y[D == 1], probs = tau_seq, na.rm = TRUE)  
Q_Y0_empirical <- quantile(Y[D == 0], probs = tau_seq, na.rm = TRUE)  
# Empirical QTE 
QTE_empirical <- Q_Y1_empirical - Q_Y0_empirical

# Plot empirical QTE
df_QTE_empirical <- data.frame(tau = tau_seq, QTE = QTE_empirical)

plot_QTE_empirical2 <- ggplot(df_QTE_empirical, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +                     
  labs(title = "",   # Empirical QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_CDF_empirical2)
print(plot_QTE_empirical2)

#------ DR-BART CDFs and QTE with CI ------ 
fit_DRBART <- drbart(Y,
                     D,
                     nburn = 2000,
                     nsim = 2000,
                     variance = "ux"
)
grid = seq(from = min(Y), to = max(Y), by = 0.05)
x_pred <- rbind(1,0)
preds_cdf <- predict(fit_DRBART, xpred = x_pred, ygrid = grid, type = "distribution")
# Extract mean CDF and 95% cred int of CDF
CDFy_drbart <- apply(preds_cdf$preds, c(1,2), mean)
CDFy_drbart_2.5 <- apply(preds_cdf$preds, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
CDFy_drbart_97.5 <- apply(preds_cdf$preds, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
# Mean CDF for Treatment and Control
CDFy1_drbart <- CDFy_drbart[1, ]
CDFy0_drbart <- CDFy_drbart[2, ]
# 2.5  CDF for Treatment and Control
CDFy1_drbart_2.5 <- CDFy_drbart_2.5[1, ]
CDFy0_drbart_2.5 <- CDFy_drbart_2.5[2, ]
# 97.5  CDF for Treatment and Control
CDFy1_drbart_97.5 <- CDFy_drbart_97.5[1, ]
CDFy0_drbart_97.5 <- CDFy_drbart_97.5[2, ]

# Plot DR-BART CDFs with cred int
group_colors <- c("Control" = "#F8766D", "Treatment" = "#00BFC4")

df_CDF_drbart <- data.frame(
  grid,                 
  cdf = c(CDFy1_drbart, CDFy0_drbart),  
  cdf_2.5 = c(CDFy1_drbart_2.5, CDFy0_drbart_2.5),
  cdf_97.5 = c(CDFy1_drbart_97.5, CDFy0_drbart_97.5),
  group = rep(c("Treatment", "Control"), each = length(grid))
)

plot_CDF_drbart2 <- ggplot(df_CDF_drbart, aes(x = grid, y = cdf, color = group, fill = group)) +
  geom_ribbon(aes(ymin = cdf_2.5, ymax = cdf_97.5), alpha = 0.2, color = NA) +  
  geom_line(size = 1) +  
  scale_color_manual(values = group_colors) +  
  scale_fill_manual(values = group_colors) +  
  labs(
    title = "",   #Estimated CDFs via DR-BART
    x = "Y",
    y = "Probability",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal()

print(plot_CDF_drbart2)

# Separate posterior samples for F1 and F0
preds1 <- preds_cdf$preds[1, , ]  
preds0 <- preds_cdf$preds[2, , ]
# Function to invert the CDF to obtain quantiles
get_quantile <- function(cdf, grid, probs) {
  sapply(probs, function(p) grid[which.min(abs(cdf - p))])  # Find closest CDF value and return corresponding Y
}
# Quantiles for each posterior sample
quant0 <- apply(preds0, 2, get_quantile, grid = grid, probs = tau_seq)  
quant1 <- apply(preds1, 2, get_quantile, grid = grid, probs = tau_seq)
# QTE for each posterior sample
qte_samples <- quant1 - quant0 
# Mean QTE 
QTE_drbart <- apply(qte_samples, 1, mean)
# Compute 95% cred int (2.5th and 97.5th percentiles)
QTE_drbart_2.5 <- apply(qte_samples, 1, quantile, probs = 0.025)
QTE_drbart_97.5 <- apply(qte_samples, 1, quantile, probs = 0.975)

# Plot DR-BART QTE with cred int
qte_drbart_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_drbart,
  QTE_2.5 = QTE_drbart_2.5,  
  QTE_97.5 = QTE_drbart_97.5 
)

plot_QTE_drbart2 <- ggplot(qte_drbart_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  #geom_hline(yintercept = ATE_drbart, linetype = "dashed", size = 1, color = "black") +  
  labs(title = "",   #DR-BART QTE and ATE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_drbart2)

#------ Firpo QTE with CI ------ 
Firpo <- ci.qte(Y ~ D,
                data=df,
                probs=tau_seq,
                se= TRUE,
                iters=100)
QTE_Firpo <- Firpo$qte
QTE_Firpo_2.5 <- Firpo$qte.lower
QTE_Firpo_97.5 <- Firpo$qte.upper 
ATE_Firpo <- Firpo$ate

# Plot Firpo QTE with conf int
qte_Firpo_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_Firpo,
  QTE_2.5 = QTE_Firpo_2.5,  
  QTE_97.5 = QTE_Firpo_97.5 
)

plot_QTE_Firpo2 <- ggplot(qte_Firpo_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  #geom_hline(yintercept = ATE_Firpo, linetype = "dashed", size = 1, color = "black") +  
  labs(title = "",   #Firpo (2007) QTE and ATE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_Firpo2)

#------ TRUE QTE and Firpo QTE with CI ------
# Plot true QTE and Firpo QTE with conf int
df_combined <- data.frame(
  tau = rep(tau_seq, 2),
  QTE = c(QTE_Firpo, QTE_true),
  Method = rep(c("Firpo (2007)", "True"), each = length(tau_seq))
)

plot_combined3 <- ggplot() +
  geom_ribbon(data = qte_Firpo_df, aes(x = tau, ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2, show.legend = FALSE) +
  geom_line(data = qte_Firpo_df, aes(x = tau, y = QTE, color = "Firpo (2007)"), size = 1) +
  geom_line(data = df_QTE_true, aes(x = Quantile, y = QTE, color = "True"), size = 0.7, linetype = "longdash") +
  labs(title = "", x = "Quantile", y = "Treatment Effect", color = "QTE") +
  scale_color_manual(values = c("True" = "red", "Firpo (2007)" = "blue")) +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_combined3)

#------ TRUE QTE and DR-BART QTE with CI ------
# Plot True QTE and DR-BART QTE with cred int
plot_combined4 <- ggplot() +
  geom_ribbon(data = qte_drbart_df, aes(x = tau, ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2, show.legend = FALSE) +
  geom_line(data = qte_drbart_df, aes(x = tau, y = QTE, color = "DR-BART"), size = 1) +
  geom_line(data = df_QTE_true, aes(x = Quantile, y = QTE, color = "True"), size = 0.7, linetype = "longdash") +
  labs(title = "", x = "Quantile", y = "Treatment Effect", color = "QTE") +
  scale_color_manual(values = c("True" = "red", "DR-BART" = "blue")) +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_combined4)

#------ RESULTS ----------
# Compute all distance measures for Empirical estimation
L1_distance_empirical <- sum(abs(QTE_true - QTE_empirical))
L2_distance_empirical <- sqrt(sum((QTE_true - QTE_empirical)^2))
Wasserstein_distance_empirical <- sum(abs(QTE_true - QTE_empirical)) / length(QTE_true)
max_abs_deviation_empirical <- max(abs(QTE_true - QTE_empirical))
ISE_empirical <- sum((QTE_true - QTE_empirical)^2) / length(QTE_true)

# Compute all distance measures for DR-BART estimation
L1_distance_drbart <- sum(abs(QTE_true - QTE_drbart))
L2_distance_drbart <- sqrt(sum((QTE_true - QTE_drbart)^2))
Wasserstein_distance_drbart <- sum(abs(QTE_true - QTE_drbart)) / length(QTE_true)
max_abs_deviation_drbart <- max(abs(QTE_true - QTE_drbart))
ISE_drbart <- sum((QTE_true - QTE_drbart)^2) / length(QTE_true)

# Compute all distance measures for Firpo estimation
L1_distance_Firpo <- sum(abs(QTE_true - QTE_Firpo))
L2_distance_Firpo <- sqrt(sum((QTE_true - QTE_Firpo)^2))
Wasserstein_distance_Firpo <- sum(abs(QTE_true - QTE_Firpo)) / length(QTE_true)
max_abs_deviation_Firpo <- max(abs(QTE_true - QTE_Firpo))
ISE_Firpo <- sum((QTE_true - QTE_Firpo)^2) / length(QTE_true)

# Results  
cat("L1 Distance DRBART:", L1_distance_drbart, "\n")
cat("L1 Distance EMPIRICAL:", L1_distance_empirical, "\n")
cat("L1 Distance Firpo:", L1_distance_Firpo, "\n")

cat("L2 Distance DRBART:", L2_distance_drbart, "\n")
cat("L2 Distance EMPIRICAL:", L2_distance_empirical, "\n")
cat("L2 Distance Firpo:", L2_distance_Firpo, "\n")

cat("Wasserstein Distance DRBART:", Wasserstein_distance_drbart, "\n")
cat("Wasserstein Distance EMPIRICAL:", Wasserstein_distance_empirical, "\n")
cat("Wasserstein Distance Firpo:", Wasserstein_distance_Firpo, "\n")

cat("Max Absolute Deviation DRBART:", max_abs_deviation_drbart, "\n")
cat("Max Absolute Deviation EMPIRICAL:", max_abs_deviation_empirical, "\n")
cat("Max Absolute Deviation Firpo:", max_abs_deviation_Firpo, "\n")

cat("ISE DRBART:", ISE_drbart, "\n\n")
cat("ISE EMPIRICAL:", ISE_empirical, "\n")
cat("ISE Firpo:", ISE_Firpo, "\n")
###################### SIMULATION 3: observational study DR-BART CQTE, QTE Gamma ###################
#------ SETUP ------- 
set.seed(42)
# X: education levels
n <- 1000
X <- sample(0:2, n, replace = TRUE, prob = c(0.3, 0.4, 0.3))  # Education levels: 0,1,2 with probs in population of 0.3,0.4,0.3
# Treatment assignment probability is based on X to mimic selection bias
p_D1_given_X <- c(0.3, 0.5, 0.7)               # Higher education -> higher chance of treatment
D <- rbinom(n, 1, prob = p_D1_given_X[X + 1])  # Treatment assignment: p(d=1|x=0)=0.3, p(d=1|x=1)=0.5, p(d=1|x=2)=0.7
# Gamma distribution parameters based on education level
shape_X <- c(2, 2.5, 3)  
scale_X <- c(1, 1.2, 1.5)  
# Generate potential outcomes
Y_0 <- rgamma(n, shape = shape_X[X + 1], scale = scale_X[X + 1])   
Y_1 <- rgamma(n, shape = shape_X[X + 1], scale = scale_X[X + 1] * 1.5)  
# Observed outcome
Y <- ifelse(D == 1, Y_1, Y_0)
# Quantile sequence
tau_seq <- seq(0.01, 0.99, by = 0.01)

#------ TRUE CQTE ------- 
df_QTE_true_X <- data.frame()
for (x in 0:2) {
  # Theoretical quantiles
  Q_Y0_true <- qgamma(tau_seq, shape = shape_X[x + 1], scale = scale_X[x + 1])
  Q_Y1_true <- qgamma(tau_seq, shape = shape_X[x + 1], scale = scale_X[x + 1] * 1.5)
  QTE_true <- Q_Y1_true - Q_Y0_true
  # True CQTE
  df_QTE_true_X <- rbind(df_QTE_true_X, data.frame(Quantile = tau_seq, QTE = QTE_true, 
                                                   Education = factor(x, levels = 0:2, labels = c("Low", "Mid", "High"))))
}

# Plot true CQTE for each education group 
plot_CQTE_true_X <- ggplot(df_QTE_true_X, aes(x = Quantile, y = QTE, color = Education)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Low" = "red3", "Mid" = "green3", "High" = "blue3")) +  
  labs(title = "",  # True CQTEs and CATEs
       x = "Quantile", y = "Treatment Effect", color = "Education Level") +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_CQTE_true_X)

#------ TRUE QTE -------
n_large <- 100000  
# X generated from the population probabilities
X_large <- sample(0:2, n_large, replace = TRUE, prob = c(0.3, 0.4, 0.3))
# Generate corresponding Y_0 and Y_1 from the Gamma distributions
Y_0_large <- rgamma(n_large, shape = shape_X[X_large + 1], scale = scale_X[X_large + 1])
Y_1_large <- rgamma(n_large, shape = shape_X[X_large + 1], scale = scale_X[X_large + 1] * 1.5)
# True QTE 
Q_Y0_true <- quantile(Y_0_large, tau_seq)  
Q_Y1_true <- quantile(Y_1_large, tau_seq)
QTE_true <- Q_Y1_true - Q_Y0_true    

# Plot true QTE
df_QTE_true <- data.frame(Quantile = tau_seq, QTE = QTE_true)

plot_QTE_true3 <- ggplot(df_QTE_true, aes(x = Quantile, y = QTE)) +
  geom_line(color = "blue", linewidth = 1) +  
  labs(title = "",  #True QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_true3)

#------ EMPIRICAL CQTE  ------- 
# Compute empirical CQTE for each education level
df_QTE_empirical_X <- data.frame()
for (x in 0:2) {
  Q_Y1_empirical <- quantile(Y[D == 1 & X == x], probs = tau_seq, na.rm = TRUE)
  Q_Y0_empirical <- quantile(Y[D == 0 & X == x], probs = tau_seq, na.rm = TRUE)
  QTE_empirical <- Q_Y1_empirical - Q_Y0_empirical
  df_QTE_empirical_X <- rbind(df_QTE_empirical_X, data.frame(Quantile = tau_seq, QTE = QTE_empirical, Education = factor(x, levels = 0:2, labels = c("Low", "Mid", "High"))))
}

# Plot Empirical CQTE for each education level
plot_CQTE_empirical_X <- ggplot(df_QTE_empirical_X, aes(x = Quantile, y = QTE, color = Education)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Low" = "red3", "Mid" = "green3", "High" = "blue3")) +  
  labs(title = "", #Empirical CQTEs and CATEs
       x = "Quantile", y = "Treatment Effect", color = "Education Level") +
  theme_minimal() + 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_CQTE_empirical_X)

#------ EMPIRICAL QTE accounting for X ------- 
# Grid of Y values for CDF estimation
Y_grid <- seq(min(Y), max(Y), length.out = 500)  
# Empirical CDFs for Y1|X and Y0|X
F_Y1_given_X <- list()
F_Y0_given_X <- list()
for (x in 0:2) {
  F_Y1_given_X[[as.character(x)]] <- ecdf(Y[D == 1 & X == x])(Y_grid)  # CDF of Y1 | X
  F_Y0_given_X[[as.character(x)]] <- ecdf(Y[D == 0 & X == x])(Y_grid)  # CDF of Y0 | X
}
# Sample probabilities P(X = x)
P_X <- table(X) / length(X)
# Weighted unconditional CDFs
F_Y1_uncond <- Reduce(`+`, mapply(`*`, F_Y1_given_X, P_X, SIMPLIFY = FALSE))  # Weighted sum
F_Y0_uncond <- Reduce(`+`, mapply(`*`, F_Y0_given_X, P_X, SIMPLIFY = FALSE))  # Weighted sum
# Function to invert the CDF to obtain quantiles
get_quantile <- function(cdf, grid, probs) {
  sapply(probs, function(p) grid[which.min(abs(cdf - p))])  # Find closest CDF value and return corresponding Y
}
# Empirical quantiles for treated and untreated groups
Q_Y1_empirical <- get_quantile(F_Y1_uncond, Y_grid, tau_seq)  # Treated
Q_Y0_empirical <- get_quantile(F_Y0_uncond, Y_grid, tau_seq)  # Untreated
QTE_empirical <- Q_Y1_empirical - Q_Y0_empirical  # Compute QTE

# Plot Empirical QTE accounting for selection bias
df_QTE_empirical <- data.frame(Quantile = tau_seq, QTE = QTE_empirical)

plot_QTE_empirical3 <- ggplot(df_QTE_empirical, aes(x = Quantile, y = QTE)) +
  geom_line(color = "blue", linewidth = 1) +  
  labs(title = "",   #Empirical QTE with CIA Adjustment
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_empirical3)

#------ DR-BART ------- 
# DR-BART on the full dataset, including treatment D and education level X as covariates
fit_DRBART <- drbart::drbart(y = Y, 
                             x = cbind(X, D),  
                             nburn = 2000, 
                             nsim = 2000, 
                             variance = "ux"  
)
# Unique values from dataset
X_levels <- sort(unique(X))  # Unique values of X (education level)
D_levels <- sort(unique(D))  # Unique values of D (treatment status)
# All possible (X, D) combinations 
x_pred <- expand.grid(X = X_levels, D = D_levels)  
x_pred <- as.matrix(x_pred)  
print(x_pred)
grid = seq(from = min(Y), to = max(Y), by = 0.05)
preds_cdf <- predict(fit_DRBART, xpred = x_pred, ygrid = grid, type = "distribution")
dim(preds_cdf$preds)

#------ DR-BART CDF -------
# Mean estimates over posterior samples
CDFy_drbart <- apply(preds_cdf$preds, c(1,2), mean)
# CDFs for each (X, D) pair
CDF_X0_D0 <- CDFy_drbart[1, ]  # Low Education, Control
CDF_X1_D0 <- CDFy_drbart[2, ]  # Mid Education, Control
CDF_X2_D0 <- CDFy_drbart[3, ]  # High Education, Control
CDF_X0_D1 <- CDFy_drbart[4, ]  # Low Education, Treated
CDF_X1_D1 <- CDFy_drbart[5, ]  # Mid Education, Treated
CDF_X2_D1 <- CDFy_drbart[6, ]  # High Education, Treated

# Plot DR-BART conditional CDFs
df_CDF <- data.frame(
  Y = rep(grid, 6),
  CDF = c(CDF_X0_D0, CDF_X0_D1, CDF_X1_D0, CDF_X1_D1, CDF_X2_D0, CDF_X2_D1),
  Education = factor(rep(c("Low", "Low", "Mid", "Mid", "High", "High"), each = length(grid)), 
                     levels = c("Low", "Mid", "High")),  # Ensure correct legend order
  Treatment = rep(c("Control", "Treatment"), times = 3, each = length(grid))
)
custom_colors <- c("Low" = "red3", "Mid" = "green3", "High" = "blue3")

plot_CDF_drbart3 <- ggplot(df_CDF, aes(x = Y, y = CDF, color = Education, linetype = Treatment)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(title = "", #Estimated CDFs via DR-BART
       x = "Y", y = "Probability", color = "Education Level", linetype = "Group") +
  theme_minimal()+ 
  theme(
    legend.position = c(0.8, 0.45),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(2,2,2,2)  
  )

print(plot_CDF_drbart3)
#------ DR-BART CQTE with CI ------ 
# Posterior samples
preds0_X0 <- preds_cdf$preds[1, , ]  
preds0_X1 <- preds_cdf$preds[2, , ]
preds0_X2 <- preds_cdf$preds[3, , ]
preds1_X0 <- preds_cdf$preds[4, , ]
preds1_X1 <- preds_cdf$preds[5, , ]
preds1_X2 <- preds_cdf$preds[6, , ]
# Function to invert the CDF to obtain quantiles
get_quantile <- function(cdf, grid, probs) {
  sapply(probs, function(p) grid[which.min(abs(cdf - p))])  # Find closest CDF value and return corresponding Y
}
# Quantiles for each posterior draw
quant0_X0 <- apply(preds0_X0, 2, get_quantile, grid = grid, probs = tau_seq)  
quant1_X0 <- apply(preds1_X0, 2, get_quantile, grid = grid, probs = tau_seq)
quant0_X1 <- apply(preds0_X1, 2, get_quantile, grid = grid, probs = tau_seq)  
quant1_X1 <- apply(preds1_X1, 2, get_quantile, grid = grid, probs = tau_seq)
quant0_X2 <- apply(preds0_X2, 2, get_quantile, grid = grid, probs = tau_seq)  
quant1_X2 <- apply(preds1_X2, 2, get_quantile, grid = grid, probs = tau_seq)
# CQTE for each posterior sample
qte_samples_X0 <- quant1_X0 - quant0_X0 
qte_samples_X1 <- quant1_X1 - quant0_X1 
qte_samples_X2 <- quant1_X2 - quant0_X2 
# Mean CQTE
QTE_drbart_X0 <- apply(qte_samples_X0, 1, mean)
QTE_drbart_X1 <- apply(qte_samples_X1, 1, mean)
QTE_drbart_X2 <- apply(qte_samples_X2, 1, mean)
# 95% cred int (2.5th and 97.5th percentiles)
QTE_drbart_X0_2.5 <- apply(qte_samples_X0, 1, quantile, probs = 0.025)
QTE_drbart_X0_97.5 <- apply(qte_samples_X0, 1, quantile, probs = 0.975)
QTE_drbart_X1_2.5 <- apply(qte_samples_X1, 1, quantile, probs = 0.025)
QTE_drbart_X1_97.5 <- apply(qte_samples_X1, 1, quantile, probs = 0.975)
QTE_drbart_X2_2.5 <- apply(qte_samples_X2, 1, quantile, probs = 0.025)
QTE_drbart_X2_97.5 <- apply(qte_samples_X2, 1, quantile, probs = 0.975)

# Plot DR-BART CQTE with cred int
cqte_est_df <- data.frame(
  tau = rep(tau_seq, 3), 
  QTE = c(QTE_drbart_X0, QTE_drbart_X1, QTE_drbart_X2), 
  QTE_2.5 = c(QTE_drbart_X0_2.5, QTE_drbart_X1_2.5, QTE_drbart_X2_2.5),  
  QTE_97.5 = c(QTE_drbart_X0_97.5, QTE_drbart_X1_97.5, QTE_drbart_X2_97.5),  
  Education = rep(c("Low", "Mid", "High"), each = length(tau_seq))  
)
cqte_est_df$Education <- factor(cqte_est_df$Education, levels = c("Low", "Mid", "High"))
education_colors <- c("Low" = "red3", "Mid" = "green3", "High" = "blue3")

plot_CQTE_drbart_CI <- ggplot(cqte_est_df, aes(x = tau, y = QTE, color = Education, fill = Education)) +
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), alpha = 0.1, color = NA) +  # Shaded 95% CI without borders
  geom_line(size = 1) +  # Plot only the mean QTE curves
  scale_color_manual(values = education_colors) +  # Set colors for QTE curves
  scale_fill_manual(values = education_colors) +  # Set fill colors for shaded CI
  labs(title = "",   #DR-BART CQTEs and CATEs
       x = "Quantile", y = "Treatment Effect", color = "Education Level", fill = "Education Level") +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_CQTE_drbart_CI)

#------ DR-BART QTE with CI ------
df <- data.frame(D=D, Y=Y, X=X)
# Empirical proportions of X education level
education_probs <- table(df$X) / nrow(df)  # P(X = x) for x = 0,1,2
p_X0 <- education_probs["0"]
p_X1 <- education_probs["1"]
p_X2 <- education_probs["2"]
# Weighted average of conditional CDFs for each simulation
preds1 <- preds1_X0 * p_X0 + preds1_X1 * p_X1 + preds1_X2 * p_X2
preds0 <- preds0_X0 * p_X0 + preds0_X1 * p_X1 + preds0_X2 * p_X2
# Function to invert the CDF to obtain quantiles
get_quantile <- function(cdf, grid, probs) {
  sapply(probs, function(p) grid[which.min(abs(cdf - p))])  # Find closest CDF value and return corresponding Y
}
# Quantiles for each posterior sample
quant0 <- apply(preds0, 2, get_quantile, grid = grid, probs = tau_seq)  
quant1 <- apply(preds1, 2, get_quantile, grid = grid, probs = tau_seq)
# QTE for each posterior sample
qte_samples <- quant1 - quant0 
# Mean QTE 
QTE_drbart <- apply(qte_samples, 1, mean)
# Compute 95% cred int (2.5th and 97.5th percentiles)
QTE_drbart_2.5 <- apply(qte_samples, 1, quantile, probs = 0.025)
QTE_drbart_97.5 <- apply(qte_samples, 1, quantile, probs = 0.975)

# Plot DR-BART QTE with CI
qte_drbart_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_drbart,
  QTE_2.5 = QTE_drbart_2.5,  
  QTE_97.5 = QTE_drbart_97.5 
)

plot_QTE_drbart3 <- ggplot(qte_drbart_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  labs(title = "",   # DR-BART QTE and ATE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_drbart3)

#------ Firpo CQTE with CI ------
df_0 <- subset(df, X == 0, select = c(D, Y))
df_1 <- subset(df, X == 1, select = c(D, Y))
df_2 <- subset(df, X == 2, select = c(D, Y))

Firpo_X0 <- ci.qte(Y ~ D,
                   data=df_0,
                   probs=tau_seq,
                   se= TRUE,
                   iters=100)
QTE_Firpo_X0 <- Firpo_X0$qte
QTE_Firpo_X0_2.5 <- Firpo_X0$qte.lower
QTE_Firpo_X0_97.5 <- Firpo_X0$qte.upper 

Firpo_X1 <- ci.qte(Y ~ D,
                   data=df_1,
                   probs=tau_seq,
                   se= TRUE,
                   iters=100)
QTE_Firpo_X1 <- Firpo_X1$qte
QTE_Firpo_X1_2.5 <- Firpo_X1$qte.lower
QTE_Firpo_X1_97.5 <- Firpo_X1$qte.upper 

Firpo_X2 <- ci.qte(Y ~ D,
                   data=df_2,
                   probs=tau_seq,
                   se= TRUE,
                   iters=100)
QTE_Firpo_X2 <- Firpo_X2$qte
QTE_Firpo_X2_2.5 <- Firpo_X2$qte.lower
QTE_Firpo_X2_97.5 <- Firpo_X2$qte.upper 

# Plot Firpo CQTE with conf int
cqte_Firpo_df <- data.frame(
  tau = rep(tau_seq, 3),  
  QTE = c(QTE_Firpo_X0, QTE_Firpo_X1, QTE_Firpo_X2),  
  QTE_2.5 = c(QTE_Firpo_X0_2.5, QTE_Firpo_X1_2.5, QTE_Firpo_X2_2.5), 
  QTE_97.5 = c(QTE_Firpo_X0_97.5, QTE_Firpo_X1_97.5, QTE_Firpo_X2_97.5),  
  Education = rep(c("Low", "Mid", "High"), each = length(tau_seq))  
)
cqte_Firpo_df$Education <- factor(cqte_Firpo_df$Education, levels = c("Low", "Mid", "High"))
education_colors <- c("Low" = "red3", "Mid" = "green3", "High" = "blue3")

plot_CQTE_Firpo_CI <- ggplot(cqte_Firpo_df, aes(x = tau, y = QTE, color = Education, fill = Education)) +
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), alpha = 0.1, color = NA) +  
  geom_line(size = 1) +  
  scale_color_manual(values = education_colors) +  
  scale_fill_manual(values = education_colors) +  
  labs(title = "",   #DR-BART CQTEs and CATEs
       x = "Quantile", y = "Treatment Effect", color = "Education Level", fill = "Education Level") +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_CQTE_Firpo_CI)
#------ Firpo QTE with CI ------ 
Firpo <- ci.qte(Y ~ D,
                xformla=~X,
                data=df,
                probs=tau_seq,
                se= TRUE,
                iters=100)
QTE_Firpo <- Firpo$qte
QTE_Firpo_2.5 <- Firpo$qte.lower
QTE_Firpo_97.5 <- Firpo$qte.upper 
ATE_Firpo <- Firpo$ate

# Plot Firpo QTE with conf int
qte_Firpo_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_Firpo,
  QTE_2.5 = QTE_Firpo_2.5,  
  QTE_97.5 = QTE_Firpo_97.5 
)

plot_QTE_Firpo3 <- ggplot(qte_Firpo_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  labs(title = "",  # Firpo (2007) QTE and ATE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_Firpo3)

#------ TRUE QTE and Firpo QTE with CI ------
# Plot TRUE QTE and Firpo QTE with conf int
df_combined <- data.frame(
  tau = rep(tau_seq, 2),
  QTE = c(QTE_Firpo, QTE_true),
  Method = rep(c("Firpo (2007)", "True"), each = length(tau_seq))
)

plot_combined5 <- ggplot() +
  geom_ribbon(data = qte_Firpo_df, aes(x = tau, ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2, show.legend = FALSE) +
  geom_line(data = qte_Firpo_df, aes(x = tau, y = QTE, color = "Firpo (2007)"), size = 1) +
  geom_line(data = df_QTE_true, aes(x = Quantile, y = QTE, color = "True"), size = 0.7, linetype = "longdash") +
  labs(title = "", x = "Quantile", y = "Treatment Effect", color = "QTE") +
  scale_color_manual(values = c("True" = "red", "Firpo (2007)" = "blue")) +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_combined5)

#------ TRUE QTE and DR-BART QTE with CI ------
# Plot True QTE and DR-BART QTE with cred int
plot_combined6 <- ggplot() +
  geom_ribbon(data = qte_drbart_df, aes(x = tau, ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2, show.legend = FALSE) +
  geom_line(data = qte_drbart_df, aes(x = tau, y = QTE, color = "DR-BART"), size = 1) +
  geom_line(data = df_QTE_true, aes(x = Quantile, y = QTE, color = "True"), size = 0.7, linetype = "longdash") +
  labs(title = "", x = "Quantile", y = "Treatment Effect", color = "QTE") +
  scale_color_manual(values = c("True" = "red", "DR-BART" = "blue")) +
  theme_minimal()+ 
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )   

print(plot_combined6)

#------ RESULTS CQTE ----------
# Extract theoretical conditional QTE
QTE_true_X0 <- df_QTE_true_X %>% filter(Education == "Low") %>% pull(QTE)
QTE_true_X1 <- df_QTE_true_X %>% filter(Education == "Mid") %>% pull(QTE)
QTE_true_X2 <- df_QTE_true_X %>% filter(Education == "High") %>% pull(QTE)

# Extract theoretical empirical conditional QTE
QTE_empirical_X0 <- df_QTE_empirical_X %>% filter(Education == "Low") %>% pull(QTE)
QTE_empirical_X1 <- df_QTE_empirical_X %>% filter(Education == "Mid") %>% pull(QTE)
QTE_empirical_X2 <- df_QTE_empirical_X %>% filter(Education == "High") %>% pull(QTE)

#### X0
# Compute all distance measures for Empirical estimation
L1_distance_empirical_X0 <- sum(abs(QTE_true_X0 - QTE_empirical_X0))
L2_distance_empirical_X0 <- sqrt(sum((QTE_true_X0 - QTE_empirical_X0)^2))
Wasserstein_distance_empirical_X0 <- sum(abs(QTE_true_X0 - QTE_empirical_X0)) / length(QTE_true_X0)
max_abs_deviation_empirical_X0 <- max(abs(QTE_true_X0 - QTE_empirical_X0))
ISE_empirical_X0 <- sum((QTE_true_X0 - QTE_empirical_X0)^2) / length(QTE_true_X0)

# Compute all distance measures for DR-BART estimation
L1_distance_drbart_X0 <- sum(abs(QTE_true_X0 - QTE_drbart_X0))
L2_distance_drbart_X0 <- sqrt(sum((QTE_true_X0 - QTE_drbart_X0)^2))
Wasserstein_distance_drbart_X0 <- sum(abs(QTE_true_X0 - QTE_drbart_X0)) / length(QTE_true_X0)
max_abs_deviation_drbart_X0 <- max(abs(QTE_true_X0 - QTE_drbart_X0))
ISE_drbart_X0 <- sum((QTE_true_X0 - QTE_drbart_X0)^2) / length(QTE_true_X0)

# Compute all distance measures for Firpo estimation
L1_distance_firpo_X0 <- sum(abs(QTE_true_X0 - QTE_Firpo_X0))
L2_distance_firpo_X0 <- sqrt(sum((QTE_true_X0 - QTE_Firpo_X0)^2))
Wasserstein_distance_firpo_X0 <- sum(abs(QTE_true_X0 - QTE_Firpo_X0)) / length(QTE_true_X0)
max_abs_deviation_firpo_X0 <- max(abs(QTE_true_X0 - QTE_Firpo_X0))
ISE_firpo_X0 <- sum((QTE_true_X0 - QTE_Firpo_X0)^2) / length(QTE_true_X0)

#### X1
# Compute all distance measures for Empirical estimation
L1_distance_empirical_X1 <- sum(abs(QTE_true_X1 - QTE_empirical_X1))
L2_distance_empirical_X1 <- sqrt(sum((QTE_true_X1 - QTE_empirical_X1)^2))
Wasserstein_distance_empirical_X1 <- sum(abs(QTE_true_X1 - QTE_empirical_X1)) / length(QTE_true_X1)
max_abs_deviation_empirical_X1 <- max(abs(QTE_true_X1 - QTE_empirical_X1))
ISE_empirical_X1 <- sum((QTE_true_X1 - QTE_empirical_X1)^2) / length(QTE_true_X1)

# Compute all distance measures for DR-BART estimation
L1_distance_drbart_X1 <- sum(abs(QTE_true_X1 - QTE_drbart_X1))
L2_distance_drbart_X1 <- sqrt(sum((QTE_true_X1 - QTE_drbart_X1)^2))
Wasserstein_distance_drbart_X1 <- sum(abs(QTE_true_X1 - QTE_drbart_X1)) / length(QTE_true_X1)
max_abs_deviation_drbart_X1 <- max(abs(QTE_true_X1 - QTE_drbart_X1))
ISE_drbart_X1 <- sum((QTE_true_X1 - QTE_drbart_X1)^2) / length(QTE_true_X1)

# Compute all distance measures for Firpo estimation
L1_distance_firpo_X1 <- sum(abs(QTE_true_X1 - QTE_Firpo_X1))
L2_distance_firpo_X1 <- sqrt(sum((QTE_true_X1 - QTE_Firpo_X1)^2))
Wasserstein_distance_firpo_X1 <- sum(abs(QTE_true_X1 - QTE_Firpo_X1)) / length(QTE_true_X1)
max_abs_deviation_firpo_X1 <- max(abs(QTE_true_X1 - QTE_Firpo_X1))
ISE_firpo_X1 <- sum((QTE_true_X1 - QTE_Firpo_X1)^2) / length(QTE_true_X1)

#### X2
# Compute all distance measures for Empirical estimation
L1_distance_empirical_X2 <- sum(abs(QTE_true_X2 - QTE_empirical_X2))
L2_distance_empirical_X2 <- sqrt(sum((QTE_true_X2 - QTE_empirical_X2)^2))
Wasserstein_distance_empirical_X2 <- sum(abs(QTE_true_X2 - QTE_empirical_X2)) / length(QTE_true_X2)
max_abs_deviation_empirical_X2 <- max(abs(QTE_true_X2 - QTE_empirical_X2))
ISE_empirical_X2 <- sum((QTE_true_X2 - QTE_empirical_X2)^2) / length(QTE_true_X2)

# Compute all distance measures for DR-BART estimation
L1_distance_drbart_X2 <- sum(abs(QTE_true_X2 - QTE_drbart_X2))
L2_distance_drbart_X2 <- sqrt(sum((QTE_true_X2 - QTE_drbart_X2)^2))
Wasserstein_distance_drbart_X2 <- sum(abs(QTE_true_X2 - QTE_drbart_X2)) / length(QTE_true_X2)
max_abs_deviation_drbart_X2 <- max(abs(QTE_true_X2 - QTE_drbart_X2))
ISE_drbart_X2 <- sum((QTE_true_X2 - QTE_drbart_X2)^2) / length(QTE_true_X2)

# Compute all distance measures for Firpo estimation
L1_distance_firpo_X2 <- sum(abs(QTE_true_X2 - QTE_Firpo_X2))
L2_distance_firpo_X2 <- sqrt(sum((QTE_true_X2 - QTE_Firpo_X2)^2))
Wasserstein_distance_firpo_X2 <- sum(abs(QTE_true_X2 - QTE_Firpo_X2)) / length(QTE_true_X2)
max_abs_deviation_firpo_X2 <- max(abs(QTE_true_X2 - QTE_Firpo_X2))
ISE_firpo_X2 <- sum((QTE_true_X2 - QTE_Firpo_X2)^2) / length(QTE_true_X2)

# Results  
cat("X0 L1 Distance DRBART:", L1_distance_drbart_X0, "\n")
cat("X0 L1 Distance EMPIRICAL:", L1_distance_empirical_X0, "\n")
cat("X0 L1 Distance FIRPO:", L1_distance_firpo_X0, "\n")
cat("X1 L1 Distance DRBART:", L1_distance_drbart_X1, "\n")
cat("X1 L1 Distance EMPIRICAL:", L1_distance_empirical_X1, "\n")
cat("X1 L1 Distance FIRPO:", L1_distance_firpo_X1, "\n")
cat("X2 L1 Distance DRBART:", L1_distance_drbart_X2, "\n")
cat("X2 L1 Distance EMPIRICAL:", L1_distance_empirical_X2, "\n")
cat("X2 L1 Distance FIRPO:", L1_distance_firpo_X2, "\n")

cat("X0 L2 Distance DRBART:", L2_distance_drbart_X0, "\n")
cat("X0 L2 Distance EMPIRICAL:", L2_distance_empirical_X0, "\n")
cat("X0 L2 Distance FIRPO:", L2_distance_firpo_X0, "\n")
cat("X1 L2 Distance DRBART:", L2_distance_drbart_X1, "\n")
cat("X1 L2 Distance EMPIRICAL:", L2_distance_empirical_X1, "\n")
cat("X1 L2 Distance FIRPO:", L2_distance_firpo_X1, "\n")
cat("X2 L2 Distance DRBART:", L2_distance_drbart_X2, "\n")
cat("X2 L2 Distance EMPIRICAL:", L2_distance_empirical_X2, "\n")
cat("X2 L2 Distance FIRPO:", L2_distance_firpo_X2, "\n")

cat("X0 Wasserstein Distance DRBART:", Wasserstein_distance_drbart_X0, "\n")
cat("X0 Wasserstein Distance EMPIRICAL:", Wasserstein_distance_empirical_X0, "\n")
cat("X0 Wasserstein Distance FIRPO:", Wasserstein_distance_firpo_X0, "\n")
cat("X1 Wasserstein Distance DRBART:", Wasserstein_distance_drbart_X1, "\n")
cat("X1 Wasserstein Distance EMPIRICAL:", Wasserstein_distance_empirical_X1, "\n")
cat("X1 Wasserstein Distance FIRPO:", Wasserstein_distance_firpo_X1, "\n")
cat("X2 Wasserstein Distance DRBART:", Wasserstein_distance_drbart_X2, "\n")
cat("X2 Wasserstein Distance EMPIRICAL:", Wasserstein_distance_empirical_X2, "\n")
cat("X2 Wasserstein Distance FIRPO:", Wasserstein_distance_firpo_X2, "\n")

cat("X0 Max Absolute Deviation DRBART:", max_abs_deviation_drbart_X0, "\n")
cat("X0 Max Absolute Deviation EMPIRICAL:", max_abs_deviation_empirical_X0, "\n")
cat("X0 Max Absolute Deviation FIRPO:", max_abs_deviation_firpo_X0, "\n")
cat("X1 Max Absolute Deviation DRBART:", max_abs_deviation_drbart_X1, "\n")
cat("X1 Max Absolute Deviation EMPIRICAL:", max_abs_deviation_empirical_X1, "\n")
cat("X1 Max Absolute Deviation FIRPO:", max_abs_deviation_firpo_X1, "\n")
cat("X2 Max Absolute Deviation DRBART:", max_abs_deviation_drbart_X2, "\n")
cat("X2 Max Absolute Deviation EMPIRICAL:", max_abs_deviation_empirical_X2, "\n")
cat("X2 Max Absolute Deviation FIRPO:", max_abs_deviation_firpo_X2, "\n")

cat("X0 ISE DRBART:", ISE_drbart_X0, "\n")
cat("X0 ISE EMPIRICAL:", ISE_empirical_X0, "\n")
cat("X0 ISE FIRPO:", ISE_firpo_X0, "\n\n")
cat("X1 ISE DRBART:", ISE_drbart_X1, "\n")
cat("X1 ISE EMPIRICAL:", ISE_empirical_X1, "\n")
cat("X1 ISE FIRPO:", ISE_firpo_X1, "\n\n")
cat("X2 ISE DRBART:", ISE_drbart_X2, "\n")
cat("X2 ISE EMPIRICAL:", ISE_empirical_X2, "\n")
cat("X2 ISE FIRPO:", ISE_firpo_X2, "\n")

#------ RESULTS QTE ----------
# Compute all distance measures between QTE_true and QTE_empirical
L1_distance_QTE_empirical <- sum(abs(QTE_true - QTE_empirical))
L2_distance_QTE_empirical <- sqrt(sum((QTE_true - QTE_empirical)^2))
Wasserstein_distance_QTE_empirical <- sum(abs(QTE_true - QTE_empirical)) / length(QTE_true)
max_abs_deviation_QTE_empirical <- max(abs(QTE_true - QTE_empirical))
ISE_QTE_empirical <- sum((QTE_true - QTE_empirical)^2) / length(QTE_true)

# Compute all distance measures between QTE_true and QTE_drbart
L1_distance_QTE_drbart <- sum(abs(QTE_true - QTE_drbart))
L2_distance_QTE_drbart <- sqrt(sum((QTE_true - QTE_drbart)^2))
Wasserstein_distance_QTE_drbart <- sum(abs(QTE_true - QTE_drbart)) / length(QTE_true)
max_abs_deviation_QTE_drbart <- max(abs(QTE_true - QTE_drbart))
ISE_QTE_drbart <- sum((QTE_true - QTE_drbart)^2) / length(QTE_true)

# Compute all distance measures between QTE_true and QTE_Firpo
L1_distance_QTE_Firpo <- sum(abs(QTE_true - QTE_Firpo))
L2_distance_QTE_Firpo <- sqrt(sum((QTE_true - QTE_Firpo)^2))
Wasserstein_distance_QTE_Firpo <- sum(abs(QTE_true - QTE_Firpo)) / length(QTE_true)
max_abs_deviation_QTE_Firpo <- max(abs(QTE_true - QTE_Firpo))
ISE_QTE_Firpo <- sum((QTE_true - QTE_Firpo)^2) / length(QTE_true)

# Results
cat("L1 Distance QTE DRBART:", L1_distance_QTE_drbart, "\n")
cat("L1 Distance QTE EMPIRICAL:", L1_distance_QTE_empirical, "\n")
cat("L1 Distance QTE FIRPO:", L1_distance_QTE_Firpo, "\n\n")

cat("L2 Distance QTE DRBART:", L2_distance_QTE_drbart, "\n")
cat("L2 Distance QTE EMPIRICAL:", L2_distance_QTE_empirical, "\n")
cat("L2 Distance QTE FIRPO:", L2_distance_QTE_Firpo, "\n\n")

cat("Wasserstein Distance QTE DRBART:", Wasserstein_distance_QTE_drbart, "\n")
cat("Wasserstein Distance QTE EMPIRICAL:", Wasserstein_distance_QTE_empirical, "\n")
cat("Wasserstein Distance QTE FIRPO:", Wasserstein_distance_QTE_Firpo, "\n\n")

cat("Max Absolute Deviation QTE DRBART:", max_abs_deviation_QTE_drbart, "\n")
cat("Max Absolute Deviation QTE EMPIRICAL:", max_abs_deviation_QTE_empirical, "\n")
cat("Max Absolute Deviation QTE FIRPO:", max_abs_deviation_QTE_Firpo, "\n\n")

cat("ISE QTE DRBART:", ISE_QTE_drbart, "\n")
cat("ISE QTE EMPIRICAL:", ISE_QTE_empirical, "\n")
cat("ISE QTE FIRPO:", ISE_QTE_Firpo, "\n")
