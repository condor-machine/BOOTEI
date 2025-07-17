# ================================================================
#   BOOTEI – Monte-Carlo study for χ² test of independence (NA-safe)
#   · two categorical variables (contingency table k × k)
#   · compare: permutation (exact) vs BOOTEI (bagged + perm)
#   · robust to n <= 10 and sparse tables
#   · repeats each replicate if any p-value is NA (up to max_attempts)
#   · supports an outer loop of n_rep independent repetitions
# ================================================================

# ────────────────────────────────────────────────────────────────
# 0.  Libraries & C++ engine
# ────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, tibble, ggplot2 …
})

library(Cairo)

Rcpp::sourceCpp("BOOTEI.cpp")
set.seed(123)

# ────────────────────────────────────────────────────────────────
# 1.  Simulation helpers (categorical × categorical)
# ────────────────────────────────────────────────────────────────

simulate_data_base <- function(n = 100, k = 2, assoc = 0.5, effe = TRUE) {
  stopifnot(k >= 2, assoc >= 0, assoc <= 1)
  
  cat_x <- LETTERS[1:k]
  cat_y <- letters[1:k]
  
  # Sample X marginally uniform
  x <- sample(cat_x, n, TRUE)
  
  if (!effe || assoc == 0) {
    y <- sample(cat_y, n, TRUE)  # independence
  } else {
    base_prob <- rep(1 / k, k)
    cond_probs <- matrix(rep(base_prob, each = k), nrow = k)
    for (i in seq_len(k)) {
      boost <- assoc * (1 - base_prob[i])
      cond_probs[i, ]  <- base_prob
      cond_probs[i, i] <- base_prob[i] + boost
      cond_probs[i, -i] <- base_prob[-i] * ((1 - cond_probs[i, i]) / sum(base_prob[-i]))
    }
    y <- character(n)
    for (i in seq_along(x)) {
      row <- match(x[i], cat_x)
      y[i] <- sample(cat_y, 1, prob = cond_probs[row, ])
    }
  }
  list(group1 = x, group2 = y)
}

simulate_data_safe <- function(..., max_retry = 100) {
  for (try in seq_len(max_retry)) {
    dat <- simulate_data_base(...)
    if (length(unique(dat$group1)) > 1 && length(unique(dat$group2)) > 1)
      return(dat)
  }
  stop("simulate_data_safe: exceeded max_retry without variability")
}

bootei_safe_p <- function(x, y, B, R) {
  tryCatch({
    bootei(x, y, test = "chisq", B = B, R = R)$p.value
  }, error = function(e) NA_real_)
}

# ────────────────────────────────────────────────────────────────
# 2.  Power function (perm vs BOOTEI)
# ────────────────────────────────────────────────────────────────

calculate_power_chisq <- function(n, k, alpha = 0.05, B = 200, R = 100, simulations = 1000,
                                  assoc = 0.5, effe = TRUE, max_attempts = 50) {
  perm_sig <- logical(simulations)
  boot_sig <- logical(simulations)
  i <- 1L
  while (i <= simulations) {
    attempts <- 0L
    repeat {
      attempts <- attempts + 1L
      if (attempts > max_attempts)
        stop("Replicate failed after ", max_attempts, " attempts (all NA)")
      
      dat <- simulate_data_safe(n = n, k = k, assoc = assoc, effe = effe)
      p_perm <- bootei_safe_p(dat$group1, dat$group2, B = 1, R = R)
      p_boot <- bootei_safe_p(dat$group1, dat$group2, B = B, R = R)
      if (anyNA(c(p_perm, p_boot))) next
      
      perm_sig[i] <- (p_perm < alpha)
      boot_sig[i] <- (p_boot < alpha)
      i <- i + 1L
      break
    }
  }
  tibble(power_perm = mean(perm_sig), power_boot = mean(boot_sig))
}

# ────────────────────────────────────────────────────────────────
# 3.  Design grid & constants
# ────────────────────────────────────────────────────────────────

n_list     <- c(5, 10, 20)
k_list     <- c(2, 5)
alpha_list <- c(0.01, 0.05, 0.10)
rho_list   <- c(0, 0.25, 0.5)

# Griglia completa
combo_grid <- expand.grid(n = n_list,
                          k = k_list,
                          alpha = alpha_list,
                          assoc = rho_list) %>%
  arrange(n, k)

# Run parameters
R_val  <- 1000
B_val  <- 100
sim_MC <- 1000  # inner MC
n_rep  <- 5    # outer reps

# ────────────────────────────────────────────────────────────────
# 4.  Main simulation loop (with progress bar)
# ────────────────────────────────────────────────────────────────

if (!requireNamespace("progress", quietly = TRUE))
  stop("Package 'progress' required for progress bar.")

total_tasks <- n_rep * nrow(combo_grid)
results <- vector("list", total_tasks)
idx <- 1L

pb <- progress::progress_bar$new(total = total_tasks, clear = FALSE,
                                 format = "[:bar] :percent | rep=:rep n=:n k=:k α=:alpha ρ=:rho | eta=:eta")

for (rep in seq_len(n_rep)) {
  for (row in seq_len(nrow(combo_grid))) {
    cfg <- combo_grid[row, ]
    n_obs <- cfg$n; k_dim <- cfg$k; a <- cfg$alpha; rho <- cfg$assoc
    pb$tick(tokens = list(rep = rep, n = n_obs, k = k_dim, alpha = a, rho = rho))
    
    res <- calculate_power_chisq(n = n_obs, k = k_dim, alpha = a, B = B_val, R = R_val,
                                 simulations = sim_MC, assoc = rho, effe = (rho > 0)) %>%
      mutate(rep = rep, n = n_obs, k = k_dim, alpha = a, assoc = rho)
    
    results[[idx]] <- res; idx <- idx + 1L
  }
}

final_df <- bind_rows(results)

# ────────────────────────────────────────────────────────────────
# 5.  Summaries & plots
# ────────────────────────────────────────────────────────────────

aggregated <- final_df %>%
  group_by(n, k, alpha, assoc) %>%
  summarise(across(starts_with("power_"), mean), .groups = "drop")

write_csv(aggregated, "C:/Users/innav/Desktop/BOOTEI/res_CHI2.csv")

aggregated <- read.csv("C:/Users/innav/Desktop/BOOTEI/res_CHI2.csv")

# 5.1 Type-I error (assoc = 0) – BAR plot

type1_plot <- aggregated %>%
  filter(assoc == 0) %>%
  pivot_longer(starts_with("power_"), names_to = "method", values_to = "emp_alpha") %>%
  mutate(alpha = factor(alpha)) %>%
  ggplot(aes(alpha, emp_alpha, fill = method)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  facet_grid(n ~ k, labeller = label_both) +
  geom_hline(aes(yintercept = as.numeric(as.character(alpha))), linetype = "dashed", colour = "grey30") +
  theme_minimal(base_size = 11) +
  labs(title = "Empirical Type-I Error – χ² (perm vs BOOTEI)",
       subtitle = "Dashed lines = nominal α", y = "Empirical α",
       x = "Nominal α", fill = "Method")

print(type1_plot)

# 5.2 Power (assoc > 0) – BAR plot stratified by α

power_plot <- aggregated %>%
  filter(assoc > 0) %>%
  mutate(assoc = factor(assoc), alpha = factor(alpha),
         k_alpha = paste0("k=", k, ", α=", alpha)) %>%
  pivot_longer(starts_with("power_"), names_to = "method", values_to = "power") %>%
  ggplot(aes(assoc, power, fill = method)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  facet_grid(n ~ k_alpha, labeller = label_both) +
  theme_minimal(base_size = 11) +
  labs(title = "Power vs Association – χ² (perm vs BOOTEI)",
       subtitle = "Columns = (k) × (α)", y = "Power",
       x = "Association ρ", fill = "Method")

print(power_plot)

# ggsave("chi2-alpha.pdf", plot = type1_plot, device = cairo_pdf, width = 6, height = 5)
# ggsave("chi2-power.pdf", plot = power_plot, device = cairo_pdf, width = 6, height = 5)
