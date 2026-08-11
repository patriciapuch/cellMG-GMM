# =============================================================================
# Simulation Study: Regularization of cellMGGMM via Condition Number Constraint
#
# Evaluates the effect of covariance regularization (maximum condition number)
# on estimation accuracy across varying sample sizes and covariance structures.
# Requires ssMRCD >= 2.0.0.
# =============================================================================

# ~10min runtime

# --- Packages and helpers ----------------------------------------------------

library(foreach)
library(doParallel)
library(doSNOW)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(gridExtra)
library(xtable)
library(tidyr)

source("sim_setup.R")  # provides genData() and eval_simul_simple()


# --- Parameter grid ----------------------------------------------------------

# Rows: (p = dimension, n = sample size per group, N = number of groups)
# p fixed at 10; n varies from under- to over-determined relative to p
pnN <- matrix(
  c(10,  8, 2,
    10, 10, 2,
    10, 12, 2,
    10, 15, 2,
    10, 20, 2,
    10, 30, 2),
  ncol     = 3,
  byrow    = TRUE,
  dimnames = list(NULL, c("p", "n", "N"))
)

cell_gamma <- 6
pi_diag    <- 1
csep       <- 0.5
cell_eps   <- c(0.1)   # clean vs 10% cellwise contamination
seed       <- 1:100
cond       <- c(100, Inf, -1)  # regularized vs unrestricted vs cellMCD
corr_type  <- c("A0X", "ALYZCOR")

pars <- expand.grid(
  cell_gamma = cell_gamma,
  pi_diag    = pi_diag,
  csep       = csep,
  cell_eps   = cell_eps,
  seed       = seed,
  corr_type  = corr_type,
  cond       = cond
)

# Combine pnN rows with the scalar parameter grid
npars <- nrow(pars)
pars <- cbind(
  matrix(rep(pnN, each = npars), byrow = FALSE, ncol = 3, nrow = npars * nrow(pnN)),
  pars
)
colnames(pars)[1:3] <- c("p", "n", "N")

message(sprintf("Total simulation runs: %d", nrow(pars)))


# --- Parallel setup ----------------------------------------------------------

n_cores <- 4
message(sprintf("Using %d parallel workers", n_cores))

cl <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(cl)

pb            <- txtProgressBar(min = 0, max = nrow(pars), style = 3)
progress      <- function(n) setTxtProgressBar(pb, n)
progress_opts <- list(progress = progress)


# --- Simulation loop ---------------------------------------------------------

results <- foreach(
  i             = 1:nrow(pars),
  .packages     = c("cellWise", "ssMRCD", "rrcov", "mclust", "tclust"),
  .options.snow = progress_opts,
  .combine      = c
) %dopar% {

  res <- tryCatch({

    simulated <- genData(
      p        = pars$p[i],
      n        = pars$n[i],
      N        = pars$N[i],
      pi_diag  = pars$pi_diag[i],
      seed     = pars$seed[i],
      cell_eps = pars$cell_eps[i],
      row_eps  = 0,                  # no rowwise outliers
      gamma    = pars$cell_gamma[i],
      cond     = 100,
      type     = pars$corr_type[i],
      csep     = pars$csep[i],
      mu_type  = "cs"
    )

    start <- Sys.time()
    if(pars$cond[i] != -1){
      out <- ssMRCD::cellMGGMM(
        X       = simulated$X,
        groups  = simulated$groups,
        alpha   = 1,
        hperc   = 0.75,
        nsteps  = 500,
        maxcond = pars$cond[i]   # 100 (regularized) or Inf (unrestricted)
      )
    }
    if(pars$cond[i] == -1 ){
      out = list()
      tmp = lapply(1:pars$N[i], function(x) cellWise::cellMCD(X = simulated$X[simulated$groups == x, ]))
      out$mu = lapply(tmp, function(x) x$mu)
      out$Sigma = lapply(tmp, function(x) x$S)
      out$W = do.call(rbind, lapply(tmp, function(x) x$W))
    }

    elapsed <- difftime(Sys.time(), start, units = "secs")

    list(out = out, simulated = simulated, time = elapsed)

  }, error = function(e) e)

  list(res)
}

close(pb)
stopCluster(cl)


# --- Save results ------------------------------------------------------------

save(results, file = "simulations_regularization_result.RData")
save(pars,    file = "simulations_regularization_parameters.RData")


# =============================================================================
# Post-processing and figure generation
# =============================================================================
load("simulations_regularization_result.RData")
load("simulations_regularization_parameters.RData")


# --- Identify failed runs ----------------------------------------------------

is_error <- sapply(results, function(x) inherits(x, "error") || "message" %in% names(x))
nonokay  <- which(is_error)
okay     <- which(!is_error)

if (length(nonokay) > 0) {
  warning(sprintf("%d run(s) failed out of %d total.", length(nonokay), nrow(pars)))
  print(apply(pars[nonokay, c("n", "cond")], MARGIN = 2, table))
  print(table(sapply(nonokay, function(x) conditionMessage(results[[x]]))))
}


# --- Labels ------------------------------------------------------------------
labels_fac <- c(
  "0"         = "mu == 0",
  "0.5"       = "mu ~ ' varying'",
  "0.9"       = "pi[diag] == 0.9",
  "0.75"      = "pi[diag] == 0.75",
  "recall"    = "'Recall'",
  "precision" = "'Precision'",
  "fscore"    = "'F1 Score'",
  "A0X"       = "A0X",
  "ALYZCOR"   = "ALYZ",
  "rho.1"     = "rho[1]",
  "rho.2"     = "rho[2]"
)

pars_scenarios <- unique(pars[okay, c("p", "N", "corr_type", "cell_eps", "pi_diag")]) %>%
  mutate(corr_type = as.character(corr_type)) %>%
  filter(corr_type == "ALYZCOR")


# --- Shared base boxplot helper ----------------------------------------------
base_boxplot <- function(dat, ylab, log_y = FALSE, ylim_fixed = NULL) {
  g <- ggplot(dat) +
    geom_boxplot(
      aes(x     = n,
          y     = value,
          group = interaction(n, cond),
          fill  = factor(cond,
                         levels = c(100, Inf, -1),
                         labels = c("Regularized", "Unregularized", "cellMCD"))),
      linewidth    = 0.3,
      median.linewidth = 0.8,
      outlier.size = 0.8,
      position = position_dodge(preserve = "single")
    ) +
    labs(x = "Sample size (n)", y = ylab, fill = "") +
    theme_bw(base_size = 16) +
    theme(legend.position = "top", panel.grid = element_blank()) +
    scale_fill_manual(values = c("lightblue4", "lightblue", "orange")) +
    scale_x_continuous(labels = c(8, 10, 12, 15, 20, 30), breaks = c(8, 10, 12, 15, 20 ,30)) +
    guides(fill = guide_legend(nrow = 1))

  if (log_y)             g <- g + scale_y_log10()
  if (!is.null(ylim_fixed)) g <- g + ylim(ylim_fixed)
  g
}


# --- Generate plots per scenario ---------------------------------------------
for (i in seq_len(nrow(pars_scenarios))) {

  scenario_key <- do.call(paste, pars_scenarios[i, 1:5, drop = FALSE])
  all_keys     <- do.call(paste, as.data.frame(pars[, c("p", "N", "corr_type", "cell_eps", "pi_diag")]))
  indices_plot <- intersect(which(all_keys == scenario_key), okay)

  if (length(indices_plot) == 0) next

  N_filter <- pars_scenarios$N[i]

  tag <- paste0(pars_scenarios$corr_type[i],
                "_eps", pars_scenarios$cell_eps[i],
                "_pi",  pars_scenarios$pi_diag[i])

  # Sigma: KL divergence
  kl_vec <- rep(NA_real_, length(indices_plot))
  for (j in seq_along(indices_plot)) {
    idx      <- indices_plot[j]
    res      <- results[[idx]]
    kl_sigma <- sapply(seq_len(N_filter),
                       function(k) KL(res$simulated$Sigma[[k]], res$out$Sigma[[k]]))
    kl_vec[j] <- mean(kl_sigma)
  }
  dat_sigma <- cbind(pars[indices_plot, ], value = kl_vec) %>%
    filter(cell_gamma %in% c(2, 6, 10), !is.na(value))

  if (nrow(dat_sigma) > 0) {
    ggsave(
      base_boxplot(dat_sigma, "KL Divergence", log_y = TRUE) + theme(aspect.ratio = 1/2),
      file = paste0("figures/simuls_reg_", tag, "_sigma.pdf"), width = 7, height = 4
    )
  }

  # mu: MSE
  mse_mu_vec <- rep(NA_real_, length(indices_plot))
  for (j in seq_along(indices_plot)) {
    idx        <- indices_plot[j]
    res        <- results[[idx]]
    mse_mu     <- sapply(seq_len(N_filter),
                         function(k) MSE_mu(res$simulated$mu[[k]], res$out$mu[[k]]))
    mse_mu_vec[j] <- mean(mse_mu)
  }
  dat_mu <- cbind(pars[indices_plot, ], value = mse_mu_vec) %>%
    filter(cell_gamma %in% c(2, 6, 10), !is.na(value))

  # probs: MSE
  mse_pi_vec <- rep(NA_real_, length(indices_plot))
  for (j in seq_along(indices_plot)) {
    idx <- indices_plot[j]
    res <- results[[idx]]
    if (any(!is.na(res$out$pi_groups))) {
      mse_pi_vec[j] <- MSE_pi(pi0   = res$simulated$pi_groups,
                              pihat = res$out$pi_groups)
    }
  }
  dat_probs <- cbind(pars[indices_plot, ], value = mse_pi_vec) %>%
    filter(cell_gamma %in% c(2, 6, 10), !is.na(value))

  if (nrow(dat_mu) > 0)    g_mu    <- base_boxplot(dat_mu,    expression(MSE(mu)), log_y = TRUE)
  if (nrow(dat_probs) > 0) g_probs <- base_boxplot(dat_probs, expression(MSE(pi))) + theme(legend.position = "none")

  if (nrow(dat_mu) > 0 && nrow(dat_probs) > 0) {
    ggsave(
      gridExtra::grid.arrange(g_mu, g_probs, ncol = 2, widths = c(2, 1)),
      file = paste0("figures/simuls_reg_", tag, "_pimu.pdf"), width = 10, height = 6
    )
  }

  # W: precision, recall, F1
  w_eval <- matrix(NA_real_, nrow = length(indices_plot),
                   ncol = 7,
                   dimnames = list(NULL, c("TPR","TNR","FPR","FNR","precision","recall","fscore")))
  for (j in seq_along(indices_plot)) {
    idx <- indices_plot[j]
    res <- results[[idx]]
    if (!is.null(res$out$W)) {
      w_eval[j, ] <- class_out(Wreal = res$simulated$W, What = res$out$W)
    }
  }
  dat_w <- cbind(pars[indices_plot, ], w_eval) %>%
    filter(cell_gamma %in% c(2, 6, 10)) %>%
    tidyr::pivot_longer(cols = TPR:fscore, names_to = "metric", values_to = "value") %>%
    filter(metric %in% c("precision", "recall", "fscore"), !is.na(value))

  if (nrow(dat_w) > 0) {
    g_w <- base_boxplot(dat_w, "", ylim_fixed = c(0, 1)) +
      facet_grid(cols     = vars(metric),
                 labeller = as_labeller(labels_fac, default = label_parsed))
    ggsave(g_w, file = paste0("figures/simuls_reg_", tag, "_W.pdf"), width = 10, height = 5)
  }

  # Replication count table
  n_reps <- dat_sigma %>%          # any dat_* with the right pars columns works
    select(n, p, N, cell_gamma, pi_diag, csep, cell_eps, seed, cond, corr_type) %>%
    distinct() %>%
    group_by(n, p, N, cell_gamma, pi_diag, csep, cell_eps, cond, corr_type) %>%
    rename(mu = csep) %>%
    mutate(mu = ifelse(mu == 0, "0", "varying"),
           cell_gamma = as.character(cell_gamma)) %>%
    summarise(N_success = n(), .groups = "drop") %>%
    as.data.frame()

  print.xtable(xtable(n_reps, caption = paste("Successful replications:", tag)),
               include.rownames = FALSE,
               file = paste0("tables/simuls_reg_", tag, "_nreps.txt"))
}


# Extract and inspect rho (applied regularization strength)

rho_mat <- do.call(rbind, lapply(results, function(x) {
  r <- x$out$rho
  if (is.null(r)) NA_real_ else r
}))
pars <- cbind(pars, rho = rho_mat)

g_rho <- ggplot(pars %>%
                  filter(cond == 100, corr_type == "ALYZCOR") %>%
                  pivot_longer(cols = c(rho.1, rho.2),
                               names_to  = "rho_index",
                               values_to = "value_rho")) +
  geom_boxplot(aes(x     = n,
                   y     = value_rho,
                   group = interaction(n)),
               linewidth    = 0.3,
               median.linewidth =0.8,
               outlier.size = 0.8) +
  facet_grid(cols     = vars(rho_index),
             labeller = as_labeller(labels_fac, default = label_parsed)) +
  labs(x = "Sample size (n)", y = expression(rho)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "top", panel.grid = element_blank())
g_rho

ggsave(g_rho, file = paste0("figures/simuls_reg_rho.pdf"), width = 10, height = 3)
