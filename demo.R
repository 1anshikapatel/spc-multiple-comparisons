# =============================================================================
# SPC Method — Quick Demo
# Example 11.5: Oil Filter Data (Math 310 textbook)
# =============================================================================

source("spc_method.R")   # loads my_comparison()

# --- Data ---
# Five brands of oil filters, J = 9 observations per group
# MSE = 0.088 from the ANOVA table
means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
J     <- 9
MSE   <- 0.088

cat("=== Oil Filter Example (Example 11.5) ===\n")
cat("Groups: 5 oil filter brands\n")
cat("Sample means:", means, "\n")
cat("J =", J, "| MSE =", MSE, "\n\n")

# --- Apply SPC method ---
results <- my_comparison(means, J, MSE, alpha = 0.05)

cat("Pairwise Comparison Results:\n")
print(results)

cat("\nSignificant pairs:", sum(results$significant), "of", nrow(results), "\n")

cat("\nNon-significant pairs (smallest differences):\n")
print(results[!results$significant, ])

cat("\nInterpretation:\n")
cat("- 8 of 10 pairs are significantly different at alpha = 0.05\n")
cat("- Groups 1 vs 4 and Groups 3 vs 5 (both diff = 0.200) are NOT significant\n")
cat("- The step-down rule stopped at rank 9; rank 10 was auto-marked non-significant\n")
