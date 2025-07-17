# ================================================================
#   BOOTEI – Monte‑Carlo study for Mann–Whitney (Wilcoxon rank‑sum)
#   · two independent ordinal samples
#   · compare: permutation (exact, BOOTEI B = 1) vs BOOTEI bagged + perm
#   · robust to n1, n2 as small as 5 and Likert scales 3–7 points
#   · repeats each replicate if any p‑value is NA (max_attempts)
#   · supports an outer loop of n_rep independent repetitions
# ================================================================

# ────────────────────────────────────────────────────────────────
# 0.  Libraries & C++ engine
# ────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
})

library(Cairo)

Rcpp::sourceCpp("BOOTEI.cpp")
set.seed(123)

# ────────────────────────────────────────────────────────────────
# 1.  Simulation helpers (independent groups)
# ────────────────────────────────────────────────────────────────

## base generator – may return constant vector
simulate_group <- function(n, k = 5, shift = 0) {
  stopifnot(n >= 2, k >= 2)
  base_prob <- rep(1 / k, k)
  if (shift == 0) {                    # uniform
    prob <- base_prob
  } else {                             # stochastically larger distribution
    w    <- (seq_len(k))^(1 + 4 * shift)  # weights grow with category index
    prob <- w / sum(w)
  }
  sample(seq_len(k), n, TRUE, prob)
}

simulate_mw_base <- function(n1 = 5, n2 = 5, k = 3, assoc = 0.5, effe = TRUE) {
  if (!effe || assoc == 0) {
    g1 <- simulate_group(n1, k, 0)
    g2 <- simulate_group(n2, k, 0)
  } else {
    g1 <- simulate_group(n1, k, 0)          # control group
    g2 <- simulate_group(n2, k, assoc)      # treatment shifted
  }
  list(group1 = g1, group2 = g2)
}

simulate_mw_safe <- function(..., max_retry = 100) {
  for (i in seq_len(max_retry)) {
    dat <- simulate_mw_base(...)
    if (length(unique(dat$group1)) > 1 && length(unique(dat$group2)) > 1)
      return(dat)
  }
  stop("simulate_mw_safe: exceeded max_retry without variability")
}

bootei_safe_p <- function(x, y, B, R) {
  tryCatch({
    bootei(x, y, test = "mannwhitney", B = B, R = R)$p.value
  }, error = function(e) NA_real_)
}

# ────────────────────────────────────────────────────────────────
# 2.  Power function (perm vs BOOTEI)
# ────────────────────────────────────────────────────────────────

calculate_power_mw <- function(n1, n2, k, alpha = 0.05, B = 200, R = 100,
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
        stop("Replicate failed after ", max_attempts, " attempts (all NA)")
      dat <- simulate_mw_safe(n1 = n1, n2 = n2, k = k, assoc = assoc, effe = effe)
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

# equal‑size scenarios
nx_equal <- c(5, 10, 20)

# unequal scenarios (n1 < n2)
extra_pairs <- tibble(n1 = c(5, 5, 10), n2 = c(10, 20, 20))

equal_pairs <- tibble(n1 = nx_equal, n2 = nx_equal)

# Combinazione completa
size_grid <- bind_rows(equal_pairs, extra_pairs) %>%
  crossing(k = c(3, 5)) 

alpha_grid <- c(0.01, 0.05, 0.10)
assoc_grid <- c(0, 0.25, 0.5)

sim_MC <- 1000   # inner MC
n_rep  <- 5     # outer reps

B_val <- 100
R_val <- 1000

# ────────────────────────────────────────────────────────────────
# 4.  Main simulation loop with progress bar
# ────────────────────────────────────────────────────────────────

if (!requireNamespace("progress", quietly = TRUE))
  stop("Package 'progress' required.")

total_tasks <- n_rep * nrow(size_grid) * length(alpha_grid) * length(assoc_grid)
results <- vector("list", total_tasks)
idx <- 1L

pb <- progress::progress_bar$new(total = total_tasks, clear = FALSE,
                                 format = "[:bar] :percent | rep=:rep n1=:n1 n2=:n2 k=:k α=:alpha ρ=:rho | eta=:eta")

for (rep in seq_len(n_rep)) {
  for (row in seq_len(nrow(size_grid))) {
    cfg <- size_grid[row, ]
    n1 <- cfg$n1; n2 <- cfg$n2; k_dim <- cfg$k
    for (a in alpha_grid) {
      for (rho in assoc_grid) {
        pb$tick(tokens = list(rep = rep, n1 = n1, n2 = n2, k = k_dim, alpha = a, rho = rho))
        res <- calculate_power_mw(n1 = n1, n2 = n2, k = k_dim, alpha = a,
                                  B = B_val, R = R_val, simulations = sim_MC,
                                  assoc = rho, effe = (rho > 0)) %>%
          mutate(rep = rep, n1 = n1, n2 = n2, k = k_dim, alpha = a, assoc = rho)
        results[[idx]] <- res; idx <- idx + 1L
      }
    }
  }
}

final_df <- bind_rows(results)

# ────────────────────────────────────────────────────────────────
# 5.  Summaries & plots
# ────────────────────────────────────────────────────────────────

aggregated <- final_df %>%
  group_by(n1, n2, k, alpha, assoc) %>%
  summarise(across(starts_with("power_"), mean), .groups = "drop")

# helper for facet label of sample sizes
aggregated <- aggregated %>%
  mutate(n_pair = paste0("n1=", n1, ", n2=", n2))

write_csv(aggregated, "C:/Users/innav/Desktop/BOOTEI/res_MW.csv")

aggregated <- read.csv("C:/Users/innav/Desktop/BOOTEI/res_MW.csv")
names(aggregated)[8] <- 'n'

# 5.1 Type‑I error (assoc = 0)

type1_plot <- aggregated %>%
  filter(assoc == 0) %>%
  pivot_longer(starts_with("power_"), names_to = "method", values_to = "emp_alpha") %>%
  mutate(alpha = factor(alpha)) %>%
  ggplot(aes(alpha, emp_alpha, fill = method)) +
  geom_col(position = position_dodge(.8), width = .7) +
  facet_grid(n ~ k, labeller = label_both) +
  geom_hline(aes(yintercept = as.numeric(as.character(alpha))), linetype = "dashed", colour = "grey30") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 6)) +
  labs(title = "Empirical Type-I Error – Man-Whitney (perm vs BOOTEI)",
       subtitle = "Dashed lines = nominal α", y = "Empirical α",
       x = "Nominal α", fill = "Method")

print(type1_plot)

# 5.2 Power (assoc > 0) – stratified by α

power_plot <- aggregated %>%
  filter(assoc > 0) %>%
  mutate(assoc = factor(assoc), alpha = factor(alpha),
         k_alpha = paste0("k=", k, ", α=", alpha)) %>%
  pivot_longer(starts_with("power_"), names_to = "method", values_to = "power") %>%
  ggplot(aes(assoc, power, fill = method)) +
  geom_col(position = position_dodge(.8), width = .7) +
  facet_grid(n ~ k_alpha, labeller = label_both) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 6)) +
  labs(title = "Power vs Association - Mann-Whitney (perm vs BOOTEI)",
       subtitle = "Rows = sample sizes, columns = (k) × (α)", y = "Power",
       x = "Association ρ", fill = "Method")

print(power_plot)

# ggsave("MW-alpha.pdf", plot = type1_plot, device = cairo_pdf, width = 6, height = 5)
# ggsave("MW-power.pdf", plot = power_plot, device = cairo_pdf, width = 6, height = 5)
