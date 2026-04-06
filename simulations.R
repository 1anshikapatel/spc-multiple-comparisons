# =============================================================================
# SPC Method — Simulations
# Part 3: FWER Validation | Part 4: Power Analysis
# =============================================================================

source("spc_method.R")   # loads my_comparison()

# =============================================================================
# PART 3 — FWER Validation
# =============================================================================
# Generate data under H0 (all means = 0), apply SPC, record whether ANY pair
# is falsely rejected. Repeat 10,000 times per scenario.
# Target: FWER between 0.03 and 0.07.

simulate_fwer <- function(I, J, sigma = 1, n_sim = 10000, alpha = 0.05) {
  false_rejections <- 0

  for (sim in 1:n_sim) {
    # Generate I groups of J observations, all from N(0, sigma^2) — H0 is true
    data <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)

    # Compute summary statistics
    means_sim <- rowMeans(data)
    MSE_sim   <- mean(apply(data, 1, var))   # pooled within-group variance

    # Apply SPC method
    res <- my_comparison(means_sim, J, MSE_sim, alpha)

    # Count simulation as a false rejection if ANY pair is declared significant
    if (any(res$significant)) {
      false_rejections <- false_rejections + 1
    }
  }

  return(false_rejections / n_sim)
}

# --- Run all 5 scenarios ---
set.seed(100)   # for reproducibility

scenarios <- data.frame(
  Scenario = c("A", "B", "C", "D", "E"),
  I        = c(3,   5,   7,   5,   5),
  J        = c(10,  10,  10,  5,   20)
)

cat("=== FWER Simulation Results (n_sim = 10,000, seed = 100) ===\n\n")

fwer_results           <- scenarios
fwer_results$k         <- choose(scenarios$I, 2)
fwer_results$n_sim     <- 10000
fwer_results$FWER_Est  <- NA
fwer_results$Target    <- "0.03 – 0.07"
fwer_results$Pass      <- NA

for (s in 1:nrow(scenarios)) {
  est <- simulate_fwer(I = scenarios$I[s], J = scenarios$J[s],
                       sigma = 1, n_sim = 10000, alpha = 0.05)
  fwer_results$FWER_Est[s] <- round(est, 4)
  fwer_results$Pass[s]     <- ifelse(est >= 0.03 & est <= 0.07, "YES", "NO")
}

print(fwer_results)
cat("\nAll scenarios pass:", all(fwer_results$Pass == "YES"), "\n")


# =============================================================================
# PART 4 — Power Analysis vs. Tukey's HSD
# =============================================================================
# Setup: I=5, J=10, sigma=1
# True means: mu1=mu2=mu3=0, mu4=mu5=delta
# Power = P(pair (1,4) declared significant)

simulate_power <- function(delta, I = 5, J = 10, sigma = 1,
                           n_sim = 1000, alpha = 0.05) {
  spc_detect   <- 0
  tukey_detect <- 0
  mu <- c(0, 0, 0, delta, delta)   # true means

  for (sim in 1:n_sim) {
    # Generate data under H1
    data <- matrix(NA, nrow = I, ncol = J)
    for (i in 1:I) data[i, ] <- rnorm(J, mean = mu[i], sd = sigma)

    means_sim <- rowMeans(data)
    MSE_sim   <- mean(apply(data, 1, var))

    # --- SPC ---
    res_spc <- my_comparison(means_sim, J, MSE_sim, alpha)
    pair_14 <- res_spc[res_spc$group1 == 1 & res_spc$group2 == 4, ]
    if (nrow(pair_14) > 0 && pair_14$significant) {
      spc_detect <- spc_detect + 1
    }

    # --- Tukey's HSD (built-in) ---
    df_long   <- data.frame(value = as.vector(t(data)),
                            group = factor(rep(1:I, each = J)))
    aov_model <- aov(value ~ group, data = df_long)
    tukey_res <- TukeyHSD(aov_model, conf.level = 1 - alpha)$group
    if ("4-1" %in% rownames(tukey_res) &&
        tukey_res["4-1", "p adj"] < alpha) {
      tukey_detect <- tukey_detect + 1
    }
  }

  list(SPC_power   = spc_detect   / n_sim,
       Tukey_power = tukey_detect / n_sim)
}

# --- Run power analysis ---
set.seed(42)
deltas <- c(0.5, 1.0, 1.5, 2.0)

cat("\n=== Power Analysis: SPC vs Tukey's HSD (I=5, J=10, n_sim=1000, seed=42) ===\n\n")

power_table <- data.frame(
  Effect_Size = deltas,
  k_pairs     = rep(choose(5, 2), length(deltas)),
  SPC_Power   = NA,
  Tukey_Power = NA,
  Difference  = NA
)

for (d in seq_along(deltas)) {
  pw <- simulate_power(deltas[d])
  power_table$SPC_Power[d]   <- round(pw$SPC_power,   4)
  power_table$Tukey_Power[d] <- round(pw$Tukey_power, 4)
  power_table$Difference[d]  <- round(pw$SPC_power - pw$Tukey_power, 3)
}

print(power_table)

# --- Power curve plot ---
plot(
  power_table$Effect_Size, power_table$SPC_Power,
  type  = "b", pch = 16, col = "blue", lwd = 2,
  ylim  = c(0, 1),
  xlab  = "Effect Size (delta)",
  ylab  = "Power (proportion of detections)",
  main  = "Power Comparison: SPC vs Tukey's HSD\n(I=5, J=10, sigma=1, 1000 simulations)"
)
lines(power_table$Effect_Size, power_table$Tukey_Power,
      type = "b", pch = 17, col = "red", lwd = 2)
abline(h = 0.8, lty = 2, col = "gray50")
text(1.7, 0.82, "0.8 benchmark", col = "gray40", cex = 0.85)
legend("topleft",
       legend = c("SPC (my method)", "Tukey's HSD"),
       col    = c("blue", "red"),
       pch    = c(16, 17), lty = 1, lwd = 2)
