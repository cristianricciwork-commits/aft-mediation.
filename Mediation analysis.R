# Simulation Parameters
set.seed(18071975)# My birth date
n <- 1000         # Sample size
a <- log(3)       # HYP triples odds of ED
b <- log(0.4)     # ED reduces median survival time by 60%
c <- log(0.5)     # HYP reduces median survival time by 50%
lambda   <- 400   # How long events take to happen
gamma    <- 1.5   # Event rate increases over the time
lambda_c <- 280   # How long censoring take to happen
gamma_c  <- 1.5   # Censoring rate increases over the time
# Data Generation
HYP <- rbinom(n, 1, 0.5)
ED_prob <- plogis(a * HYP)
ED <- rbinom(n, 1, ED_prob)
T <- lambda * exp(c * HYP + b * ED) * (-log(runif(n)))^(1 / gamma)
C <- rweibull(n, shape = gamma_c, scale = lambda_c)
time <- pmin(T, C)
status <- as.numeric(T <= C)
# Data Merging
df <- data.frame(HYP, ED, time, status)

# Fit AFT models
library(survival)
aft_WB <- survreg(Surv(time, status) ~ HYP + ED, data = df, dist = "weibull")
aft_LN <- survreg(Surv(time, status) ~ HYP + ED, data = df, dist = "lognormal")
aft_LL <- survreg(Surv(time, status) ~ HYP + ED, data = df, dist = "loglogistic")
# Compare model fits
AIC(aft_WB, aft_LN, aft_LL)
BIC(aft_WB, aft_LN, aft_LL)

#ESTIMATION OF THE TIME TO EVENT
# Fit AFT WB model
aft_WB <- survreg(Surv(time, status) ~ HYP + ED, data = df, dist = "weibull")
# Extract coefficients
coeffs <- coef(aft_WB)
scale <- aft_WB$scale  

#Function to compute median survival time
M_time <- function(HYP, ED, coeffs, scale) {
    mu <- coeffs["(Intercept)"] + coeffs["HYP"] * HYP + coeffs["ED"] * ED
    median_time <- exp(mu) * (log(2))^scale
    return(median_time)}

#Compute for all combinations
combinations <- expand.grid(HYP = c(0, 1), ED = c(0, 1))
combinations$MedianTime <- apply(combinations, 1, function(row) {
M_time(HYP = row["HYP"], ED = row["ED"], coeffs, scale)})

#Print the results
print(combinations)

#AFT MEDIATION ANALYSIS CUSTOM PROGRAM WITH BOOTSTRAP
#Part I: Settings
set.seed(21082023)  # Allegra's birth date
library(survival)
n_boot <- 5000
n <- nrow(df)
boot_results <- matrix(NA, nrow = n_boot, ncol = 3)
colnames(boot_results) <- c("NDE", "NIE", "TE")
#PART II: Bootstrap loop
for (i in 1:n_boot) {
#PART II-A: Modelling
  idx <- sample(1:n, size = n, replace = TRUE)
  df_boot <- df[idx, ]
  # Mediator model
  med_fit <- tryCatch(
    glm(ED ~ HYP, family = binomial(), data = df_boot),
    error = function(e) NULL)
  # Outcome model
  aft_fit <- tryCatch(
    survreg(Surv(time, status) ~ HYP + ED, data = df_boot, dist = "weibull"),
    error = function(e) NULL)
  # If either model failed, store NA and continue
  if (is.null(med_fit) || is.null(aft_fit)) {
    boot_results[i, ] <- NA
    next}
#Part II-B: Extraction of coefficients and scale
  coefs <- coef(aft_fit)
  scale <- aft_fit$scale
  # Median prediction function
  predict_median <- function(HYP_val, ED_val) {
    mu <- coefs["(Intercept)"] + coefs["HYP"] * HYP_val + coefs["ED"] * ED_val
    exp(mu) * (log(2))^scale}
  # Mediator probabilities
  p_ED_H0 <- predict(med_fit, newdata = data.frame(HYP = 0), type = "response")
  p_ED_H1 <- predict(med_fit, newdata = data.frame(HYP = 1), type = "response")
  # Counter factual medians
  T00 <- predict_median(0, 0)
  T01 <- predict_median(0, 1)
  T10 <- predict_median(1, 0)
  T11 <- predict_median(1, 1)
  Y0 <- (1 - p_ED_H0) * T00 + p_ED_H0 * T01
  Y1_star <- (1 - p_ED_H0) * T10 + p_ED_H0 * T11
  Y1 <- (1 - p_ED_H1) * T10 + p_ED_H1 * T11
#PART II-C: Estimation of Counter factual effects
  NDE <- Y1_star - Y0
  NIE <- Y1 - Y1_star
  TE <- Y1 - Y0
  boot_results[i, ] <- c(NDE, NIE, TE)}
# Calculate 95% CIs ignoring NA bootstrap samples
CI_lower <- apply(boot_results, 2, quantile, probs = 0.025, na.rm = TRUE)
CI_upper <- apply(boot_results, 2, quantile, probs = 0.975, na.rm = TRUE)
# Calculate point estimates on full sample (original data)
model.m <- glm(ED ~ HYP, family = binomial(), data = df)
model.y <- survreg(Surv(time, status) ~ HYP + ED, data = df, dist = "weibull")
coefs <- coef(model.y)
scale <- model.y$scale
predict_median <- function(HYP_val, ED_val) {
  mu <- coefs["(Intercept)"] + coefs["HYP"] * HYP_val + coefs["ED"] * ED_val
  exp(mu) * (log(2))^scale
}
#PART III: Output of results
p_ED_H0 <- predict(model.m, newdata = data.frame(HYP = 0), type = "response")
p_ED_H1 <- predict(model.m, newdata = data.frame(HYP = 1), type = "response")
T00 <- predict_median(0, 0)
T01 <- predict_median(0, 1)
T10 <- predict_median(1, 0)
T11 <- predict_median(1, 1)
Y0 <- (1 - p_ED_H0) * T00 + p_ED_H0 * T01
Y1_star <- (1 - p_ED_H0) * T10 + p_ED_H0 * T11
Y1 <- (1 - p_ED_H1) * T10 + p_ED_H1 * T11
NDE <- Y1_star - Y0
NIE <- Y1 - Y1_star
TE <- Y1 - Y0

# Print results with bootstrap CIs
cat("Counterfactual Mediation Analysis with 95% Bootstrap CIs:\n")
cat(sprintf("Baseline (H=0):        %.1f\n", Y0))
cat(sprintf("Natural Direct Effect: %.1f (%.1f, %.1f)\n", NDE, CI_lower["NDE"], CI_upper["NDE"]))
cat(sprintf("Natural Indirect Effect: %.1f (%.1f, %.1f)\n", NIE, CI_lower["NIE"], CI_upper["NIE"]))
cat(sprintf("Total Effect:          %.1f (%.1f, %.1f)\n", TE, CI_lower["TE"], CI_upper["TE"]))



