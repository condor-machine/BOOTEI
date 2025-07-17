# ================================================================
#   BOOTEI – Monte‑Carlo study for Spearman correlation (NA‑safe)
#   · two ordinal (Likert) variables
#   · compare: permutation (exact) vs BOOTEI (bagged + perm)
#   · robust to n <= 10 and discrete Likert scales
#   · repeats each replicate if any p‑value is NA (up to max_attempts)
#   · supports an outer loop of n_rep independent repetitions
# ================================================================

# ────────────────────────────────────────────────────────────────
# 0.  Libraries & C++ engine
# ────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)        # mvrnorm for latent normal simulation
  library(Cairo)
})

Rcpp::sourceCpp("C:/Users/innav/Desktop/BOOTEI/BOOTEI.cpp")
set.seed(123)

# ────────────────────────────────────────────────────────────────
# 1.  Safe helpers – simulation & p‑values (no parametric)
# ────────────────────────────────────────────────────────────────

simulate_likert_base <- function(n = 100, kx = 5, ky = 5, assoc = 0.5, effe = TRUE) {
  stopifnot(n >= 2, kx >= 2, ky >= 2, assoc >= 0, assoc <= 1)
  if (!effe || assoc == 0) {
    x <- sample(kx, n, TRUE)
    y <- sample(ky, n, TRUE)
  } else {
    Sigma <- matrix(c(1, assoc, assoc, 1), 2)
    z     <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
    cuts_x <- qnorm(seq(0, 1, length.out = kx + 1))
    cuts_y <- qnorm(seq(0, 1, length.out = ky + 1))
    x <- cut(z[, 1], breaks = cuts_x, labels = FALSE, include.lowest = TRUE)
    y <- cut(z[, 2], breaks = cuts_y, labels = FALSE, include.lowest = TRUE)
  }
  list(group1 = x, group2 = y)
}

simulate_likert_safe <- function(..., max_retry = 100) {
  for (i in seq_len(max_retry)) {
    dat <- simulate_likert_base(...)
    if (length(unique(dat$group1)) > 1 && length(unique(dat$group2)) > 1)
      return(dat)
  }
  stop("simulate_likert_safe: exceeded max_retry without variability")
}

bootei_safe_p <- function(x, y, B, R) {
  tryCatch({
    bootei(x, y, test = "spearman", B = B, R = R)$p.value
  }, error = function(e) NA_real_)
}

# ────────────────────────────────────────────────────────────────
# 2.  Power function (perm vs BOOTEI)
# ────────────────────────────────────────────────────────────────

calculate_power_spear <- function(n, kx, ky, alpha = 0.05, B = 200, R = 100,
                                  simulations = 1000, assoc = 0.5, effe = TRUE,
                                  max_attempts = 50) {
  perm_sig <- logical(simulations)
  boot_sig <- logical(simulations)
  i <- 1L
  while (i <= simulations) {
    attempts <- 0L
    repeat {
      attempts <- attempts + 1L
      if (attempts > max_attempts)
        stop("Replicate failed after ", max_attempts,
             " attempts (all NA)")
      
      dat <- simulate_likert_safe(n = n, kx = kx, ky = ky, assoc = assoc, effe = effe)
      p_perm <- bootei_safe_p(dat$group1, dat$group2, B = 1,   R = R)
      p_boot <- bootei_safe_p(dat$group1, dat$group2, B = B,   R = R)
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

n <- c(5, 10, 20)
k <- c(3, 5)

design_grid <- expand.grid(n = n, k = k) %>%
  arrange(n, k)

alpha_grid <- c(0.01, 0.05, 0.10)
assoc_grid <- c(0, 0.25, 0.5)

sim_MC <- 1000   # inner Monte‑Carlo size
n_rep  <- 5     # outer repetitions

B_val <- 100
R_val <- 1000

# ────────────────────────────────────────────────────────────────
# 4.  Main Monte‑Carlo loop
# ────────────────────────────────────────────────────────────────

if (!requireNamespace("progress", quietly = TRUE))
  stop("Package 'progress' required.")

total_tasks <- n_rep * nrow(design_grid) * length(alpha_grid) * length(assoc_grid)
results <- vector("list", total_tasks)
idx <- 1L

pb <- progress::progress_bar$new(total = total_tasks, clear = FALSE,
                                 format = "[:bar] :percent | rep=:rep n=:n k=:k α=:alpha ρ=:rho | eta=:eta")

for (rep in seq_len(n_rep)) {
  for (dg in seq_len(nrow(design_grid))) {
    n_obs <- design_grid$n[dg]
    k_dim <- design_grid$k[dg]
    for (a in alpha_grid) {
      for (rho in assoc_grid) {
        pb$tick(tokens = list(rep = rep, n = n_obs, k = k_dim, alpha = a, rho = rho))
        res <- calculate_power_spear(n = n_obs, kx = k_dim, ky = k_dim,
                                     alpha = a, B = B_val, R = R_val,
                                     simulations = sim_MC, assoc = rho,
                                     effe = (rho > 0)) %>%
          mutate(rep = rep, n = n_obs, k = k_dim, alpha = a, assoc = rho)
        results[[idx]] <- res
        idx <- idx + 1L
      }
    }
  }
}

final_df <- bind_rows(results)

# ────────────────────────────────────────────────────────────────
# 5.  Summaries & plots (perm vs BOOTEI)
# ────────────────────────────────────────────────────────────────

aggregated <- final_df %>%
  group_by(n, k, alpha, assoc) %>%
  summarise(across(starts_with("power_"), mean), .groups = "drop")

write_csv(aggregated, "C:/Users/innav/Desktop/BOOTEI/res_SPEARMAN.csv")

aggregated <- read.csv("C:/Users/innav/Desktop/BOOTEI/res_SPEARMAN.csv")

# 5.1 Type‑I error (assoc = 0)

type1_tbl <- aggregated %>% filter(assoc == 0)

# Bar plot

type1_plot <- type1_tbl %>%
  pivot_longer(starts_with("power_"), names_to = "method", values_to = "emp_alpha") %>%
  mutate(alpha = factor(alpha)) %>%
  ggplot(aes(alpha, emp_alpha, fill = method)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  facet_grid(n ~ k, labeller = label_both) +
  geom_hline(aes(yintercept = as.numeric(as.character(alpha))),
             linetype = "dashed", colour = "grey30") +
  theme_minimal(base_size = 11) +
  labs(title = "Empirical Type-I Error (perm vs BOOTEI)",
       subtitle = "Dashed lines = nominal α", y = "Empirical α",
       x = "Nominal α", fill = "Method")

print(type1_plot)

# 5.2 Power (assoc > 0) – stratified by α

power_tbl <- aggregated %>% filter(assoc > 0)

power_plot <- power_tbl %>%
  mutate(assoc  = factor(assoc),
         alpha  = factor(alpha),
         k_alpha = paste0("k=", k, ", α=", alpha)) %>%
  pivot_longer(starts_with("power_"), names_to = "method", values_to = "power") %>%
  ggplot(aes(assoc, power, fill = method)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  facet_grid(n ~ k_alpha, labeller = label_both) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 6)) +
  labs(title = "Power vs Association – perm vs BOOTEI",
       subtitle = "Columns = (k) × (α)", y = "Power",
       x = "Association ρ", fill = "Method")

print(power_plot)

# ggsave("C:/Users/innav/Desktop/BOOTEI/spearman-alpha.pdf", plot = type1_plot, device = cairo_pdf, width = 6, height = 5)
# ggsave("C:/Users/innav/Desktop/BOOTEI/spearman-power.pdf", plot = power_plot, device = cairo_pdf, width = 6, height = 5)
