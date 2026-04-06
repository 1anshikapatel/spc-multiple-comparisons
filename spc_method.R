# =============================================================================
# Method:      Sequential Power-Correction (SPC)
# Author:      [Your Name]
# Date:        April 2026
# Description: A step-down multiple comparison procedure for post-hoc ANOVA
#              testing. Assigns each pairwise comparison an adjusted alpha
#              based on rank, using exponent 4/5 to balance FWER control
#              (~0.05) with higher power than Tukey's HSD. Testing stops at
#              the first non-significant rank.
# =============================================================================

my_comparison <- function(means, J, MSE, alpha = 0.05) {
  # --------------------------------------------------------------------------
  # Inputs:
  #   means : numeric vector of I sample means (one per group)
  #   J     : number of observations per group (balanced design)
  #   MSE   : mean square error from the ANOVA table
  #   alpha : desired family-wise error rate (default 0.05)
  #
  # Output:
  #   data frame with columns:
  #     group1, group2  — indices of the two groups being compared
  #     diff            — absolute difference |x̄_i − x̄_j|
  #     significant     — TRUE if the pair is declared significantly different
  # --------------------------------------------------------------------------

  # --- Step 1: Basic ANOVA quantities --------------------------------------
  I        <- length(means)      # number of groups
  k        <- choose(I, 2)       # total pairwise comparisons: I(I-1)/2
  df_error <- I * (J - 1)        # error degrees of freedom from ANOVA
  se_diff  <- sqrt(2 * MSE / J)  # standard error of any pairwise difference
                                  # (estimated by plugging MSE in for sigma^2)

  # --- Step 2: Enumerate all pairs and compute absolute differences ---------
  group1 <- integer(k)
  group2 <- integer(k)
  diffs  <- numeric(k)

  idx <- 1
  for (i in 1:(I - 1)) {
    for (j in (i + 1):I) {
      group1[idx] <- i
      group2[idx] <- j
      diffs[idx]  <- abs(means[i] - means[j])
      idx <- idx + 1
    }
  }

  pairs <- data.frame(
    group1      = group1,
    group2      = group2,
    diff        = diffs,
    significant = rep(FALSE, k)   # default: not significant
  )

  # --- Step 3: Rank pairs from largest to smallest difference ---------------
  order_idx    <- order(pairs$diff, decreasing = TRUE)
  pairs_ranked <- pairs[order_idx, ]

  # --- Step 4: Apply SPC decision rule at each rank r ----------------------
  # Adjusted alpha for rank r:  alpha_r = alpha / (k - r + 1)^(4/5)
  #
  #   r = 1 (largest diff)  → (k)^(4/5) in denominator → STRICTEST threshold
  #   r = k (smallest diff) → (1)^(4/5) = 1 in denominator → most LENIENT
  #
  # Exponent 4/5 is less than 1 (Bonferroni/Holm) but greater than 0 (no
  # correction), giving a middle-ground correction that improves power while
  # keeping FWER within [0.03, 0.07].
  #
  # Step-down rule: if rank r fails, STOP — all remaining ranks stay FALSE.

  for (r in 1:k) {
    alpha_r       <- alpha / ((k - r + 1)^(4/5))  # adjusted alpha for rank r
    t_crit        <- qt(1 - alpha_r / 2, df_error) # two-sided critical t-value
    critical_diff <- t_crit * se_diff               # minimum detectable difference

    if (pairs_ranked$diff[r] > critical_diff) {
      pairs_ranked$significant[r] <- TRUE    # significant → continue to next rank
    } else {
      break                                  # not significant → stop all remaining
    }
  }

  # --- Step 5: Return results sorted by group indices ----------------------
  result <- pairs_ranked[order(pairs_ranked$group1, pairs_ranked$group2), ]
  rownames(result) <- NULL
  return(result[, c("group1", "group2", "diff", "significant")])
}
