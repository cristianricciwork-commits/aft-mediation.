# Simulation Parameters
set.seed(18071975)  # My birth date
n <- 1000           # Sample size
a <- log(3)         # HYP triples odds of ED
b <- log(0.4)       # ED reduces median survival time by 60%
c <- log(0.5)       # HYP reduces median survival time by 50%
A <- 400            # Explicit shape parameter 
gamma <- 1.5        # Shape parameter for CVD events
lambda_c <- 280     # Scale parameter for censor
gamma_c <- 1.5      # Shape parameter for censor
# Data generation
HYP <- rbinom(n, 1, 0.5)    #Subjects with HYP
ED_prob <- plogis(a * HYP)
ED <- rbinom(n, 1, ED_prob) #Subjects with ED
T <- A * exp(c * HYP + b * ED) * (-log(runif(n)))^(1 / gamma)
C <-rweibull(n,shape = gamma_c, scale=lambda_c)
time <-pmin(T,C) 
status <- as.numeric(T<=C)
# Data Merging
df <- data.frame(HYP, ED, time, status)

#Fitting distribution
library(survival)
AFT_WB <-survreg(Surv(time,status) ~ HYP + ED, data=df, dist = "weibull")
AFT_LN <-survreg(Surv(time,status) ~ HYP + ED, data=df, dist = "lognorma")
AFT_LL <-survreg(Surv(time,status) ~ HYP + ED, data=df, dist = "loglogistic")

AIC (AFT_WB,AFT_LN,AFT_LL)
BIC (AFT_WB,AFT_LN,AFT_LL)

#CUSTOM MEDIATION ANALYSIS for AFT
library(survival)

set.seed(12345)
n_boot <- 2000

calc_effects <- function(data) {
  aft_fit <- survreg(Surv(time, status) ~ HYP + ED, data = data, dist = "weibull")
  beta <- coef(aft_fit)
  scale <- aft_fit$scale
  
  predict_median_time <- function(HYP_val, ED_val) {
    lp <- beta["(Intercept)"] + beta["HYP"] * HYP_val + beta["ED"] * ED_val
    exp(lp) * (log(2))^scale
  }
  
  p_ED_given_HYP <- tapply(data$ED, data$HYP, mean)
  n <- nrow(data)
  
  median_T_01 <- mean(
    sapply(1:n, function(i) {
      p1 <- p_ED_given_HYP["1"]
      ED_sim <- rbinom(1, 1, p1)
      predict_median_time(0, ED_sim)
    })
  )
  
  median_T_10 <- mean(
    sapply(1:n, function(i) {
      p0 <- p_ED_given_HYP["0"]
      ED_sim <- rbinom(1, 1, p0)
      predict_median_time(1, ED_sim)
    })
  )
  
  median_T_00 <- mean(
    sapply(1:n, function(i) {
      p0 <- p_ED_given_HYP["0"]
      ED_sim <- rbinom(1, 1, p0)
      predict_median_time(0, ED_sim)
    })
  )
  
  median_T_11 <- mean(
    sapply(1:n, function(i) {
      p1 <- p_ED_given_HYP["1"]
      ED_sim <- rbinom(1, 1, p1)
      predict_median_time(1, ED_sim)
    })
  )
  
  NIE <- median_T_01 - median_T_00
  NDE <- median_T_10 - median_T_00
  TE  <- median_T_11 - median_T_00
  
  return(c(NIE = NIE, NDE = NDE, TE = TE, T0 = median_T_00))
}

# Point estimates
point_estimates <- calc_effects(df)

# Bootstrap results matrix (NIE, NDE, TE, T0)
boot_results <- replicate(n_boot, {
  idx <- sample(seq_len(nrow(df)), replace = TRUE)
  boot_df <- df[idx, ]
  tryCatch(calc_effects(boot_df), error = function(e) rep(NA, 4))
})

# Remove bootstrap samples with NAs
valid_idx <- complete.cases(t(boot_results))
boot_results <- boot_results[, valid_idx]

# Extract bootstrap results separately
boot_NIE <- boot_results["NIE", ]
boot_NDE <- boot_results["NDE", ]
boot_TE  <- boot_results["TE", ]
boot_T0  <- boot_results["T0", ]

# Calculate 95% CIs for raw effects
CI_lower <- apply(boot_results[1:3, ], 1, quantile, probs = 0.025)
CI_upper <- apply(boot_results[1:3, ], 1, quantile, probs = 0.975)

# Calculate % reductions (negative effect means reduction)
pct_NIE <- 100 * point_estimates["NIE"] / point_estimates["T0"]
pct_NDE <- 100 * point_estimates["NDE"] / point_estimates["T0"]
pct_TE  <- 100 * point_estimates["TE"]  / point_estimates["T0"]

# Calculate bootstrap CIs for % reductions
boot_pct_NIE <- 100 * boot_NIE / boot_T0
boot_pct_NDE <- 100 * boot_NDE / boot_T0
boot_pct_TE  <- 100 * boot_TE  / boot_T0

CI_pct_lower <- c(
  quantile(boot_pct_NIE, 0.025),
  quantile(boot_pct_NDE, 0.025),
  quantile(boot_pct_TE,  0.025)
)

CI_pct_upper <- c(
  quantile(boot_pct_NIE, 0.975),
  quantile(boot_pct_NDE, 0.975),
  quantile(boot_pct_TE,  0.975)
)

names(CI_pct_lower) <- names(CI_pct_upper) <- c("NIE", "NDE", "TE")

# Print results
cat("Estimated effects with 95% bootstrap CIs:\n")
cat(sprintf("Baseline median survival (T0): %.3f units\n\n", point_estimates["T0"]))

cat(sprintf("NIE: %.3f (%.3f, %.3f) units (%.1f%% reduction, %.1f%% to %.1f%% CI)\n",
            point_estimates["NIE"], CI_lower["NIE"], CI_upper["NIE"],
            abs(pct_NIE), abs(CI_pct_lower["NIE"]), abs(CI_pct_upper["NIE"])))

cat(sprintf("NDE: %.3f (%.3f, %.3f) units (%.1f%% reduction, %.1f%% to %.1f%% CI)\n",
            point_estimates["NDE"], CI_lower["NDE"], CI_upper["NDE"],
            abs(pct_NDE), abs(CI_pct_lower["NDE"]), abs(CI_pct_upper["NDE"])))

cat(sprintf("TE : %.3f (%.3f, %.3f) units (%.1f%% reduction, %.1f%% to %.1f%% CI)\n",
            point_estimates["TE"], CI_lower["TE"], CI_upper["TE"],
            abs(pct_TE), abs(CI_pct_lower["TE"]), abs(CI_pct_upper["TE"])))

#COUNTERFACTUAL MEDIATION ANALYSIS with mediation (Gold Standard)
library(mediation)
library(survival)

# Mediator model (logistic regression)
model.m <- glm(ED ~ HYP, family = binomial(), data = df)

# Outcome model (Weibull AFT model)
model.y <- survreg(Surv(time, status) ~ HYP + ED, data = df, dist = "weibull")

# Mediation analysis with outcome specified
results.med <- mediate(model.m,
                       model.y,
                       treat = "HYP",
                       mediator = "ED",
                       outcome = "time",  
                       boot = TRUE,
                       sims = 5000)
# Summary of results
summary(results.med)
