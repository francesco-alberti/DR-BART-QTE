library(drbart)     
library(ggplot2)     
library(qte)         
library(tidyverse)  
library(dplyr)       
library(reshape2)         
library(weights)

###################### APPLICATION 1: randomized experiment DR-BART QTE Lalonde (1986) ######################  
#------ Dataset -------
# Experimental part of Lalonde dataset
dataexp <- lalonde.exp
# Experimental data: outcome y and treatment d
Y <- dataexp$re78
D <- dataexp$treat
summary(Y)
table(D)
#------ DR-BART ------- 
fit_DRBART <- drbart::drbart(Y, 
                             D,
                             nburn = 1000,
                             nsim = 1000,
                             nthin = 1,
                             variance = "ux"
)
grid = seq(from = min(Y), to = max(Y), by = 60)
x_pred <- rbind(1,0)
preds_cdf_DRBART <- predict(fit_DRBART, xpred = x_pred, ygrid = grid, type = "distribution")
preds_cdf <- preds_cdf_DRBART
#------ DR-BART CDF with CI ------- 
# Mean CDF and 95% cred int
CDF_drbart <- apply(preds_cdf$preds, c(1,2), mean)
CDF_drbart_2.5 <- apply(preds_cdf$preds, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
CDF_drbart_97.5 <- apply(preds_cdf$preds, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
# Mean CDF for Treatment and Control
CDF1_drbart <- CDF_drbart[1, ]
CDF0_drbart <- CDF_drbart[2, ]
# 2.5  CDF for Treatment and Control
CDF1_drbart_2.5 <- CDF_drbart_2.5[1, ]
CDF0_drbart_2.5 <- CDF_drbart_2.5[2, ]
# 97.5  CDF for Treatment and Control
CDF1_drbart_97.5 <- CDF_drbart_97.5[1, ]
CDF0_drbart_97.5 <- CDF_drbart_97.5[2, ]

# Plot DR-BART CDF with cred int
group_colors <- c("Control" = "#F8766D", "Treatment" = "#00BFC4")

df_CDF_drbart <- data.frame(
  grid,                 
  cdf = c(CDF1_drbart, CDF0_drbart),  
  cdf_2.5 = c(CDF1_drbart_2.5, CDF0_drbart_2.5),
  cdf_97.5 = c(CDF1_drbart_97.5, CDF0_drbart_97.5),
  group = rep(c("Treatment", "Control"), each = length(grid))
)

plot_CDF_drbart4 <- ggplot(df_CDF_drbart, aes(x = grid, y = cdf, color = group, fill = group)) +
  geom_ribbon(aes(ymin = cdf_2.5, ymax = cdf_97.5), alpha = 0.2, color = NA) +  
  geom_line(size = 1) +  
  scale_color_manual(values = group_colors) +  
  scale_fill_manual(values = group_colors) +  
  labs(
    title = "",  #Estimated CDFs via DR-BART
    x = "Y",
    y = "Probability",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal()+ 
  theme(
    legend.position = c(0.8, 0.2),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )
print(plot_CDF_drbart4)
#------ DR-BART QTE with CI ------ 
# Separate posterior samples for Y_1 and Y_0
preds1 <- preds_cdf$preds[1, , ]  
preds0 <- preds_cdf$preds[2, , ]
# Quantile sequence
tau_seq <- seq(0.05, 0.95, by = 0.01)
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
# 95% cred int (2.5th and 97.5th percentiles)
QTE_drbart_2.5 <- apply(qte_samples, 1, quantile, probs = 0.025)
QTE_drbart_97.5 <- apply(qte_samples, 1, quantile, probs = 0.975)

# Plot DR-BART QTE with cred int
qte_est_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_drbart,
  QTE_2.5 = QTE_drbart_2.5,  
  QTE_97.5 = QTE_drbart_97.5 
)

plot_QTE_drbart4 <- ggplot(qte_est_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  labs(title = "", #DR-BART QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_drbart4)

#------ Firpo QTE with CI ------- 
df <- data.frame(D = D, Y = Y)
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
qte_est_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_Firpo,
  QTE_2.5 = QTE_Firpo_2.5,  
  QTE_97.5 = QTE_Firpo_97.5 
)

plot_QTE_Firpo4 <- ggplot(qte_est_df, aes(x = tau, y = QTE)) +
  geom_line(color = "red", size = 1) +  
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "red", alpha = 0.2) +  
  labs(title = "", #Firpo (2007) QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_Firpo4)

#------ DR-BART QTE with CI and Firpo QTE with CI ------
# Plot DR-BART QTE with cred int and Firpo QTE with conf int
combined_qte_ci_df <- data.frame(
  tau = rep(tau_seq, 2),
  QTE = c(QTE_drbart, QTE_Firpo),
  QTE_2.5 = c(QTE_drbart_2.5, QTE_Firpo_2.5),
  QTE_97.5 = c(QTE_drbart_97.5, QTE_Firpo_97.5),
  Method = rep(c("DR-BART", "Firpo (2007)"), each = length(tau_seq))
)

plot_combined_ci_dashed <- ggplot() +
  geom_line(data = combined_qte_ci_df, aes(x = tau, y = QTE, color = Method), size = 1) +
  geom_line(data = subset(combined_qte_ci_df, Method == "Firpo (2007)"), 
            aes(x = tau, y = QTE_2.5), color = "red", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qte_ci_df, Method == "Firpo (2007)"), 
            aes(x = tau, y = QTE_97.5), color = "red", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qte_ci_df, Method == "DR-BART"), 
            aes(x = tau, y = QTE_2.5), color = "blue", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qte_ci_df, Method == "DR-BART"), 
            aes(x = tau, y = QTE_97.5), color = "blue", linetype = "dashed", size = 0.25) +
  scale_color_manual(values = c("DR-BART" = "blue", "Firpo (2007)" = "red")) +
  labs(title = "",  # QTE: DR-BART vs Firpo (2007) with CI
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal() + 
  theme(
    legend.position = c(0.15, 0.8),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )

print(plot_combined_ci_dashed)

###################### APPLICATION 2: observational study DR-BART QTE Lalonde (1986) ######################  
#------ Dataset ------- 
# Observational part of Lalonde dataset
datapsid <- lalonde.psid
# Observational data: outcome y, treatment d, age, education
Y <- datapsid$re78
D <- datapsid$treat
Age <- datapsid$age
Edu <- datapsid$education
df <- data.frame(D=D, Y=Y, Edu=Edu, Age=Age)
# Create Age categories with numeric values: 
# 0 1st quartile, 1 between 1st and 3rd quartile, 2 above 3rd quartile
df$A <- cut(df$Age, 
            breaks = c(-Inf, 25, 43.5, Inf), 
            labels = c(0, 1, 2), 
            right = TRUE)
df$A <- as.numeric(as.character(df$A))
# Create Education categories with numeric values: 
# 0 1st quartile, 1 between 1st and 3rd quartile, 2 above 3rd quartile
df$E <- cut(df$Edu, 
            breaks = c(-Inf, 10, 14, Inf), 
            labels = c(0, 1, 2), 
            right = TRUE)
df$E <- as.numeric(as.character(df$E))
#------ DR-BART ------- 
# DR-BART on the full dataset, including education level E,age level A, treatment D 
fit_DRBART2 <- drbart::drbart(y = df$Y, 
                              x = cbind(df$E,df$A,df$D),  
                              nburn = 1000, 
                              nsim = 1000, 
                              nthin = 1, 
                              variance = "ux"  
)
# Extract unique values from dataset
edu_levels <- sort(unique(df$E))  
age_levels <- sort(unique(df$A))
treatment_levels <- sort(unique(df$D))  
# All possible (E,A,D) combinations 
x_pred <- expand.grid(E = edu_levels, A = age_levels,D = treatment_levels)  
x_pred <- as.matrix(x_pred)
print(x_pred)
grid = seq(from = min(Y), to = max(Y), by = 120)
preds_cdf <- predict(fit_DRBART2, xpred = x_pred, ygrid = grid, type = "distribution")
###################### QTE ###################### 
#------ DR-BART QTE with CI -------
# Extract conditional CDFs manually
preds0_E0_A0 <- preds_cdf$preds[1, , ] 
preds0_E1_A0 <- preds_cdf$preds[2, , ] 
preds0_E2_A0 <- preds_cdf$preds[3, , ] 
preds0_E0_A1 <- preds_cdf$preds[4, , ] 
preds0_E1_A1 <- preds_cdf$preds[5, , ] 
preds0_E2_A1 <- preds_cdf$preds[6, , ] 
preds0_E0_A2 <- preds_cdf$preds[7, , ] 
preds0_E1_A2 <- preds_cdf$preds[8, , ] 
preds0_E2_A2 <- preds_cdf$preds[9, , ] 

preds1_E0_A0 <- preds_cdf$preds[10, , ] 
preds1_E1_A0 <- preds_cdf$preds[11, , ] 
preds1_E2_A0 <- preds_cdf$preds[12, , ] 
preds1_E0_A1 <- preds_cdf$preds[13, , ] 
preds1_E1_A1 <- preds_cdf$preds[14, , ] 
preds1_E2_A1 <- preds_cdf$preds[15, , ] 
preds1_E0_A2 <- preds_cdf$preds[16, , ] 
preds1_E1_A2 <- preds_cdf$preds[17, , ] 
preds1_E2_A2 <- preds_cdf$preds[18, , ] 

# Compute Empirical Joint Probabilities P(A, E) 
joint_probs <- prop.table(table(df$A, df$E))
# Conditional CDFs for control group D=0
preds0_list <- list(
  "0_0" = preds0_E0_A0, "1_0" = preds0_E1_A0, "2_0" = preds0_E2_A0,
  "0_1" = preds0_E0_A1, "1_1" = preds0_E1_A1, "2_1" = preds0_E2_A1,
  "0_2" = preds0_E0_A2, "1_2" = preds0_E1_A2, "2_2" = preds0_E2_A2
)
# Conditional CDFs for treatment group D=1
preds1_list <- list(
  "0_0" = preds1_E0_A0, "1_0" = preds1_E1_A0, "2_0" = preds1_E2_A0,
  "0_1" = preds1_E0_A1, "1_1" = preds1_E1_A1, "2_1" = preds1_E2_A1,
  "0_2" = preds1_E0_A2, "1_2" = preds1_E1_A2, "2_2" = preds1_E2_A2
)
preds0 <- matrix(0, nrow = nrow(preds0_E0_A0), ncol = ncol(preds0_E0_A0))
preds1 <- matrix(0, nrow = nrow(preds1_E0_A0), ncol = ncol(preds1_E0_A0))
# Compute Weighted Average to get preds0 and preds1 
for (a in age_levels) {
  for (e in edu_levels) {
    key <- paste(a, e, sep = "_")  
    prob <- joint_probs[a + 1, e + 1]  
    if (!is.na(prob) && key %in% names(preds0_list) && key %in% names(preds1_list)) {
      preds0 <- preds0 + preds0_list[[key]] * prob  
      preds1 <- preds1 + preds1_list[[key]] * prob  
    }
  }
}
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
qte_est_df <- data.frame(
  tau = tau_seq,
  QTE = QTE_drbart,
  QTE_2.5 = QTE_drbart_2.5,  
  QTE_97.5 = QTE_drbart_97.5 
)

plot_QTE_drbart5 <- ggplot(qte_est_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  labs(title = "",  #DR-BART QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()
print(plot_QTE_drbart5)

#------ Firpo QTE with CI -------
Firpo <- ci.qte(Y ~ D,
                xformla= ~ A + E,
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

plot_QTE_Firpo5 <- ggplot(qte_Firpo_df, aes(x = tau, y = QTE)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTE_2.5, ymax = QTE_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  labs(title = "",  #Firpo (2007) QTE
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTE_Firpo5)

#------ DR-BART QTE with CI and Firpo QTE with CI ------
# Plot DR-BART QTE with cred int and Firpo QTE with conf int
combined_qte_ci_df <- data.frame(
  tau = rep(tau_seq, 2),
  QTE = c(QTE_drbart, QTE_Firpo),
  QTE_2.5 = c(QTE_drbart_2.5, QTE_Firpo_2.5),
  QTE_97.5 = c(QTE_drbart_97.5, QTE_Firpo_97.5),
  Method = rep(c("DR-BART", "Firpo (2007)"), each = length(tau_seq))
)

plot_combined_ci_dashed2 <- ggplot() +
  geom_line(data = combined_qte_ci_df, aes(x = tau, y = QTE, color = Method), size = 1) +
  geom_line(data = subset(combined_qte_ci_df, Method == "Firpo (2007)"), 
            aes(x = tau, y = QTE_2.5), color = "red", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qte_ci_df, Method == "Firpo (2007)"), 
            aes(x = tau, y = QTE_97.5), color = "red", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qte_ci_df, Method == "DR-BART"), 
            aes(x = tau, y = QTE_2.5), color = "blue", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qte_ci_df, Method == "DR-BART"), 
            aes(x = tau, y = QTE_97.5), color = "blue", linetype = "dashed", size = 0.25) +
  # Define color scheme for methods
  scale_color_manual(values = c("DR-BART" = "blue", "Firpo (2007)" = "red")) +
  labs(title = "",  #QTE: DR-BART vs Firpo (2007) with CI
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal() + 
  theme(
    legend.position = c(0.8, 0.8),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )

print(plot_combined_ci_dashed2)
###################### QTET ###################### 
#------ DR-BART QTET with CI -------
# Empirical Joint Probabilities P(A, E | D=1)
joint_probs_treated <- prop.table(table(df$A[df$D == 1], df$E[df$D == 1]))  
# Conditional CDFs for control group (D=0, counterfactual for treated)
preds0_list <- list(
  "0_0" = preds0_E0_A0, "1_0" = preds0_E1_A0, "2_0" = preds0_E2_A0,
  "0_1" = preds0_E0_A1, "1_1" = preds0_E1_A1, "2_1" = preds0_E2_A1,
  "0_2" = preds0_E0_A2, "1_2" = preds0_E1_A2, "2_2" = preds0_E2_A2
)
# Conditional CDFs for treatment group (D=1, actual treated outcomes)
preds1_list <- list(
  "0_0" = preds1_E0_A0, "1_0" = preds1_E1_A0, "2_0" = preds1_E2_A0,
  "0_1" = preds1_E0_A1, "1_1" = preds1_E1_A1, "2_1" = preds1_E2_A1,
  "0_2" = preds1_E0_A2, "1_2" = preds1_E1_A2, "2_2" = preds1_E2_A2
)
# Initialize Empty Matrices for unconditional CDFs 
# The shape is same as each conditional CDF
preds0_treated <- matrix(0, nrow = nrow(preds0_E0_A0), ncol = ncol(preds0_E0_A0))  # 1010 × 500
preds1_treated <- matrix(0, nrow = nrow(preds1_E0_A0), ncol = ncol(preds1_E0_A0))  # 1010 × 500
# Compute Weighted Average to Get preds0_treated and preds1_treated
for (a in age_levels) {
  for (e in edu_levels) {
    key <- paste(a, e, sep = "_") 
    prob <- joint_probs_treated[a + 1, e + 1] 
        if (!is.na(prob) && key %in% names(preds0_list) && key %in% names(preds1_list)) {
      preds0_treated <- preds0_treated + preds0_list[[key]] * prob  
      preds1_treated <- preds1_treated + preds1_list[[key]] * prob 
    }
  }
}
# Function to invert the CDF to obtain quantiles
get_quantile <- function(cdf, grid, probs) {
  sapply(probs, function(p) grid[which.min(abs(cdf - p))])  
}
# Quantiles for each posterior draw
quant0_treated <- apply(preds0_treated, 2, get_quantile, grid = grid, probs = tau_seq)  
quant1_treated <- apply(preds1_treated, 2, get_quantile, grid = grid, probs = tau_seq)
# QTET for each posterior sample
qtet_samples <- quant1_treated - quant0_treated 
# Mean QTET
QTET_drbart <- apply(qtet_samples, 1, mean)
# 95% credible interval (2.5th and 97.5th percentiles)
QTET_drbart_2.5 <- apply(qtet_samples, 1, quantile, probs = 0.025, na.rm = TRUE)
QTET_drbart_97.5 <- apply(qtet_samples, 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot DR-BART QTET wich cred int
qtet_est_df <- data.frame(
  tau = tau_seq,
  QTET = QTET_drbart,
  QTET_2.5 = QTET_drbart_2.5,  
  QTET_97.5 = QTET_drbart_97.5 
)

plot_QTET_drbart <- ggplot(qtet_est_df, aes(x = tau, y = QTET)) +
  geom_line(color = "blue", size = 1) +  # Estimated QTE curve
  geom_ribbon(aes(ymin = QTET_2.5, ymax = QTET_97.5), fill = "blue", alpha = 0.2) +  # Shaded 95% CI region
  labs(title = "",   #DR-BART QTET
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTET_drbart)
#------ Firpo QTET with CI -------
Firpo2 <- ci.qtet(Y ~ D,
                  xformla= ~ A + E,
                  data=df,
                  probs=tau_seq,
                  se= TRUE,
                  iters=100)
QTET_Firpo <- Firpo2$qte
QTET_Firpo_2.5 <- Firpo2$qte.lower
QTET_Firpo_97.5 <- Firpo2$qte.upper 

# Plot Firpo QTET with conf int
qtet_Firpo_df <- data.frame(
  tau = tau_seq,
  QTET = QTET_Firpo,
  QTET_2.5 = QTET_Firpo_2.5,  
  QTET_97.5 = QTET_Firpo_97.5 
)

plot_QTET_Firpo <- ggplot(qtet_Firpo_df, aes(x = tau, y = QTET)) +
  geom_line(color = "blue", size = 1) +  
  geom_ribbon(aes(ymin = QTET_2.5, ymax = QTET_97.5), fill = "blue", alpha = 0.2) + 
  labs(title = "", #Firpo (2007) QTET
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal()

print(plot_QTET_Firpo)
#------ DR-BART QTET with CI and Firpo QTET with CI ------
# Plot DR-BART QTET with cred int and Firpo QTET with conf int
combined_qtet_ci_df <- data.frame(
  tau = rep(tau_seq, 2),
  QTET = c(QTET_drbart, QTET_Firpo),
  QTET_2.5 = c(QTET_drbart_2.5, QTET_Firpo_2.5),
  QTET_97.5 = c(QTET_drbart_97.5, QTET_Firpo_97.5),
  Method = rep(c("DR-BART", "Firpo (2007)"), each = length(tau_seq))
)

plot_combined_qtet_ci_dashed <- ggplot() +
  geom_line(data = combined_qtet_ci_df, aes(x = tau, y = QTET, color = Method), size = 1) +
  geom_line(data = subset(combined_qtet_ci_df, Method == "Firpo (2007)"), 
            aes(x = tau, y = QTET_2.5), color = "red", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qtet_ci_df, Method == "Firpo (2007)"), 
            aes(x = tau, y = QTET_97.5), color = "red", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qtet_ci_df, Method == "DR-BART"), 
            aes(x = tau, y = QTET_2.5), color = "blue", linetype = "dashed", size = 0.25) +
  geom_line(data = subset(combined_qtet_ci_df, Method == "DR-BART"), 
            aes(x = tau, y = QTET_97.5), color = "blue", linetype = "dashed", size = 0.25) +
  scale_color_manual(values = c("DR-BART" = "blue", "Firpo (2007)" = "red")) +
  labs(title = "",  #QTET: DR-BART vs Firpo (2007) with CI
       x = "Quantile", y = "Treatment Effect") +
  theme_minimal() + 
  theme(
    legend.position = c(0.8, 0.8),  
    legend.background = element_rect(fill = "white", color = "white"), 
    legend.margin = margin(5,5,5,5)  
  )

print(plot_combined_qtet_ci_dashed)
