# =============================================================================
# Simulation Study: Hyperparameter tuning and sensitivity
# =============================================================================


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
library(ggh4x)

source("sim_setup.R")


# --- Parameter grid ----------------------------------------------------------
pnN <- matrix(
  c(10,  100, 2),
  ncol     = 3,
  byrow    = TRUE,
  dimnames = list(NULL, c("p", "n", "N"))
)

cell_gamma <- c(0, 2, 6, 10)
pi_diag    <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
csep       <- c(0,0.5)
cell_eps   <- c(0.1)
seed       <- 1:100
corr_type  <- c("A0X", "ALYZCOR")
alpha      <- seq(0.5, 1, 0.05)

pars <- expand.grid(
  cell_gamma = cell_gamma,
  pi_diag    = pi_diag,
  csep       = csep,
  cell_eps   = cell_eps,
  seed       = seed,
  corr_type  = corr_type,
  alpha       = alpha
)

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
  .packages     = c("cellWise", "ssMRCD"),
  .options.snow = progress_opts,
  .combine      = c
) %dopar% {

  res <- tryCatch({

    # set contamination
    if(pars$cell_gamma[i] == 0) {
      cell_eps = 0
    } else {
      cell_eps = 0.1
    }

    simulated <- genData(
      p        = pars$p[i],
      n        = pars$n[i],
      N        = pars$N[i],
      pi_diag  = pars$pi_diag[i],
      seed     = pars$seed[i],
      cell_eps = cell_eps,
      row_eps  = 0,                  # no rowwise outliers
      gamma    = pars$cell_gamma[i],
      cond     = 100,
      type     = pars$corr_type[i],
      csep     = pars$csep[i],
      mu_type  = "cs"
    )

    start <- Sys.time()
    out <- ssMRCD::cellMGGMM(
      X       = simulated$X,
      groups  = simulated$groups,
      alpha   = pars$alpha[i],
      hperc   = 0.75,
      nsteps  = 500,
      maxcond = 100
    )

    elapsed <- difftime(Sys.time(), start, units = "secs")
    list(out = out, simulated = simulated, time = elapsed)

  }, error = function(e) e)

  list(res)
}

close(pb)
stopCluster(cl)


# --- Save results ------------------------------------------------------------

save(results, file = "simulations_alpha_result.RData")
save(pars,    file = "simulations_alpha_parameters.RData")


# =============================================================================
# Post-processing and figure generation
# =============================================================================

load("simulations_alpha_result.RData")
load("simulations_alpha_parameters.RData")

# Add scenario column
pars = mutate(pars,
              scenario = case_when(p == 10 & n == 100 & N == 2 ~ "balanced2",
                                   p == 10 & n == -2 ~ "unbalanced",
                                   p == 60 ~ "highdim",
                                   p == 20 ~ "mediumdim",
                                   p == 10 & n == 100 & N == 5 ~ "balanced5"),
              scenario_nr = case_when(scenario == "balanced2" ~ "Scenario 1",
                                      scenario == "balanced5" ~ "Scenario 2",
                                      scenario == "unbalanced" ~ "Scenario 3",
                                      scenario == "mediumdim" ~ "Scenario 4",
                                      scenario == "highdim" ~ "Scenario 5",
                                      TRUE ~ NA_character_),
              corr_type_name = case_when(corr_type == "A0X" ~ "Toeplitz",
                                         corr_type == "ALYZCOR" ~ "Agostinelli et al. (2015)"),
              cell_gamma_char =  factor(paste0("cd", cell_gamma), ordered = TRUE, levels = c("cd0", "cd2", "cd4", "cd6", "cd8", "cd10")))

labels_facets = c( "0" = "mu == 0",
                   "0.5" = "mu ~ ' varying'",
                   "0.9" = "pi[diag] == 0.9",
                   "0.75" = "pi[diag] == 0.75",
                   "1" = "pi[diag] == 1",
                   "cd0" = "'No contamination'",
                   "cd2" = "gamma[cell] == 2",
                   "cd4" = "gamma[cell] == 4",
                   "cd6" = "'Contamination'",
                   "cd8" = "gamma[cell] == 8",
                   "cd10" = "gamma[cell] == 10",
                   "recall" = "'Recall'",
                   "precision" = "'Precision'",
                   "fscore" = "'F1 Score'",
                   "MSE_all" = "'All'",
                   "MSE_noise" = "'Noisy labels'",
                   "MSE_clean" = "'Clean labels'",
                   "ARI_class" = "'ARI'",
                   "ARI_diff" = "'ARI (difference)'",
                   "ARI_group" = "'ARI (pre-defined)'")


#### TEST FOR FAILED RUNS ####
# Check for non-successful runs, dependence on parameters and error messages
nonokay = which(sapply(results, function(x) "message" %in% names(x)))
okay = which(!1:dim(pars)[1] %in% nonokay)
if(length(nonokay) != 0) warning("Not all runs were successful!")

#### SENSITIVITY ANALYSIS ####
basic_plot_tile = function(data, title, trafo = "log", accuracy = 1){

  g = ggplot(data = data) +
    geom_tile(aes(y = alpha,
                  x = pi_diag,
                  fill = mean,
                  col = mean,
                  group = interaction(p, n, N, pi_diag, csep, cell_eps, corr_type, cell_gamma_char)),
              alpha = 1) +
    facet_grid(rows = vars(csep),
               cols = vars(cell_gamma_char),
               labeller = as_labeller(labels_facets,
                                      default = label_parsed),
               scales = "free_x"
    ) +
    geom_line(data = data.frame(x = c(0.475, 1.025), y = c(0.475, 1.025)),
              aes(x = x, y = y),
              col = "grey") +
    labs(title = title,
         y = expression(alpha),
         x = expression(pi[diag])
    ) +
    scale_color_distiller(trans = trafo, na.value = "white", palette = "Greens", labels = scales::number_format(accuracy = accuracy)) +
    scale_fill_distiller(trans = trafo, na.value = "white", palette = "Greens", labels = scales::number_format(accuracy = accuracy)) +
    theme_bw(base_size = 16) +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          aspect.ratio = 1/2
    ) #+
    # guides(
    #   colour = guide_colorbar(barwidth = 15),
    #   fill   = guide_colorbar(barwidth = 15)
    # )
  return(g)
}

combinations = unique(pars[, c("corr_type", "scenario")])
for(i in 1:nrow(combinations)){

  index_subset = which(pars$corr_type == combinations[i, 1] & pars$scenario == combinations[i, 2] & pars$cell_gamma %in% c(0, 6))
  N_filter = 2
  if(combinations[i, 2] == "balanced5") N_filter = 5


  #### Covariance KL ####
  kl_plot = rep(NA, length(index_subset))
  for(j in 1:length(index_subset)){
    jj = index_subset[j]
    kl_sigma = sapply(1:length(results[[jj]]$simulated$Sigma), FUN = function(x) KL(results[[jj]]$simulated$Sigma[[x]], results[[jj]]$out$Sigma[[x]]))
    kl_plot[j] = mean(kl_sigma)
  }
  dat_plot = mutate(pars[index_subset, ], value = kl_plot) %>%
    group_by(p, n, N, cell_gamma_char, pi_diag, csep, cell_eps, corr_type, alpha) %>%
    summarise(mean = mean(value), sd = sd(value))
  g = basic_plot_tile(dat_plot, "", "log") +
    labs(colour = "KL-Divergence\n", fill = "KL-Divergence\n")
  g

  ggsave(g, file = paste0("figures/simuls_alpha_", combinations$scenario[i], "_", combinations$corr_type[i], "_sigma.pdf"),
         width = 9, height = 5)



  #### Mean and pi ####
  mse_mu_plot = rep(NA, length(index_subset))
  mse_pi_plot = rep(NA, length(index_subset))

  for(j in 1:length(index_subset)){
    jj = index_subset[j]
    mse_mu = sapply(1:length(results[[jj]]$simulated$Sigma),
                    FUN = function(x) MSE_mu(results[[jj]]$simulated$mu[[x]], results[[jj]]$out$mu[[x]]))
    mse_mu_plot[j] = mean(mse_mu)

    if(any(!is.na(results[[jj]]$out$pi_groups))){
      mse_pi_plot[j] =  MSE_pi(pi0 = results[[jj]]$simulated$pi_groups,
                               pihat = results[[jj]]$out$pi_groups)
    } else mse_pi_plot[j] = NA
  }

  dat_plot = mutate(pars[index_subset, ], value = mse_mu_plot) %>%
    group_by(p, n, N, cell_gamma_char, pi_diag, csep, cell_eps, corr_type, alpha) %>%
    summarise(mean = mean(value), sd = sd(value))
  g_mu = basic_plot_tile(dat_plot, "", "log", 0.01)  +
    labs(colour = expression(atop(MSE(mu), "")), fill = expression(atop(MSE(mu), "")))
  g_mu

  ggsave(g_mu, file = paste0("figures/simuls_alpha_", combinations$scenario[i], "_", combinations$corr_type[i], "_mu.pdf"), width = 9, height = 5)

  dat_plot = mutate(pars[index_subset, ], value = mse_pi_plot) %>%
    group_by(p, n, N, cell_gamma_char, pi_diag, csep, cell_eps, corr_type, alpha) %>%
    summarise(mean = mean(value), sd = sd(value))
  g_pi = basic_plot_tile(dat_plot, "", "identity", 0.01)  +
    labs(colour = expression(atop(MSE(pi), "")), fill = expression(atop(MSE(pi), "")))
  g_pi

  ggsave(g_pi, file = paste0("figures/simuls_alpha_", combinations$scenario[i], "_", combinations$corr_type[i], "_pi.pdf"), width = 9, height = 5)


  #### W ####
  w_eval = matrix(NA,
                  ncol = 7,
                  nrow = length(index_subset),
                  dimnames = list(NULL, c("TPR","TNR","FPR", "FNR","precision", "recall","fscore")))
  for(j in 1:length(index_subset)){
    jj = index_subset[j]
    w_eval[j, ] = if(!is.null(results[[jj]]$out$W)){
      class_out(Wreal = results[[jj]]$simulated$W,
                What = results[[jj]]$out$W)
    } else {rep(NA, 7)}
  }
  dat_plot = cbind(pars[index_subset, ], w_eval) %>%
    tidyr::pivot_longer(cols = TPR:fscore) %>%
    filter(name %in% c("precision", "recall","fscore")) %>%
    group_by(p, n, N, cell_gamma_char, pi_diag, csep, cell_eps, corr_type, alpha, name) %>%
    summarise(mean = mean(value), sd = sd(value))

  g_w = basic_plot_tile(dat_plot, "", "identity", 0.1) +
    ggh4x::facet_nested(name + csep ~  cell_gamma_char,
                 labeller = as_labeller(labels_facets, default = label_parsed)) +
    labs(colour = "", fill = "")
  g_w

  ggsave(g_w, file = paste0("figures/simuls_alpha_", combinations$scenario[i], "_", combinations$corr_type[i], "_W.pdf"), width = 9, height = 9)




  #### Class Probabilities ####
  t_eval = matrix(NA,
                  ncol = 6,
                  nrow = length(index_subset),
                  dimnames = list(NULL, c("ARI_group", "ARI_class", "ARI_diff", "MSE_all", "MSE_noise", "MSE_clean")))
  for(j in 1:length(index_subset)){
    jj = index_subset[j]
    if(!is.null(results[[jj]]$out$probs) & any(!is.na(results[[jj]]$out$probs))){
      class = apply(X = results[[jj]]$out$probs, FUN = which.max, MARGIN = 1)
      t_eval[j, 1:3] = ARI(source = results[[jj]]$simulated$source,
                           groups = results[[jj]]$simulated$groups,
                           class = class)
      t_eval[j, 4:6] = MSE_probs(source = results[[jj]]$simulated$source,
                                 that = results[[jj]]$out$probs,
                                 groups = results[[jj]]$simulated$groups)[1:3]
    } else {
      t_eval[j, ] = rep(NA, 6)
    }
  }

  dat_plot = cbind(pars[index_subset, ], t_eval) %>%
    tidyr::pivot_longer(cols = ARI_group:MSE_clean) %>%
    filter(name %in% c("ARI_class")) %>%
    group_by(p, n, N, cell_gamma_char, pi_diag, csep, cell_eps, corr_type, alpha, name) %>%
    summarise(mean = mean(value), sd = sd(value))
  g_that_ari = basic_plot_tile(dat_plot, "", "identity", 0.1) +
    # facet_nested(name + csep ~ cell_gamma_char,
    #              labeller = as_labeller(labels_facets,
    #                                     default = label_parsed)) +
    labs(colour = "ARI", fill = "ARI")
  g_that_ari

  ggsave(g_that_ari, file = paste0("figures/simuls_alpha_", combinations$scenario[i], "_", combinations$corr_type[i], "_ARIW.pdf"), width = 9, height = 5)



  dat_plot = cbind(pars[index_subset, ], t_eval) %>%
    tidyr::pivot_longer(cols = ARI_group:MSE_clean) %>%
    filter(name %in% c("MSE_all", "MSE_noise")) %>%
    group_by(p, n, N, cell_gamma_char, pi_diag, csep, cell_eps, corr_type, alpha, name) %>%
    summarise(mean = mean(value), sd = sd(value))
  g_that_mse = basic_plot_tile(dat_plot, "", "identity", 0.1) +
    facet_nested(name + csep ~ cell_gamma_char,
                 labeller = as_labeller(labels_facets,
                                        default = label_parsed)) +
    labs(colour = "MSE\n", fill = "MSE\n")
  g_that_mse

  ggsave(g_that_mse, file = paste0("figures/simuls_alpha_", combinations$scenario[i], "_", combinations$corr_type[i], "_MSEW.pdf"), width = 9, height = 11)



  #### HYPERPARAMETER TUNING ####
  values_unique = unique(pars[index_subset, !colnames(pars) %in% c("alpha", "time")])
  alphastar = rep(NA, nrow(values_unique))
  for(j in 1:nrow(values_unique)){

    index_alpha = which(pars$p == values_unique$p[j] &
                        pars$n == values_unique$n[j] &
                        pars$N == values_unique$N[j] &
                        pars$cell_gamma == values_unique$cell_gamma[j] &
                        pars$pi_diag == values_unique$pi_diag[j] &
                        pars$csep == values_unique$csep[j] &
                        pars$cell_eps == values_unique$cell_eps[j] &
                        pars$seed == values_unique$seed[j] &
                        pars$corr_type == values_unique$corr_type[j] &
                        pars$scenario == values_unique$scenario[j] )
    alphastar[j] = select_alpha_simulations(list_results = results[index_alpha], alpha_grid = pars$alpha[index_alpha])[1]
    if(is.na(alphastar[j])) stop()
  }

  dat_plot = mutate(values_unique, value = alphastar)
  g_alpha = ggplot(data = dat_plot) +
    geom_boxplot(aes(y = value, x = pi_diag,
                     group = pi_diag),
               ) +
    geom_line(data = data.frame(x = c(0.5, 1), y = c(0.5, 1)),
              aes(x = x, y = y), col = "grey") +
    facet_grid(rows = vars(csep),
               cols = vars(cell_gamma_char),
               labeller = as_labeller(labels_facets,
                                      default = label_parsed),
               scales = "free_x"
    ) +
    labs(title = "",
         y = expression(alpha),
         x = expression(pi[diag])
    ) +
    theme_bw(base_size = 16) +
    theme(legend.position = "top",
          panel.grid = element_blank(),
          aspect.ratio = 1

    )
  g_alpha
  ggsave(g_alpha, file = paste0("figures/simuls_alpha_", combinations$scenario[i], "_", combinations$corr_type[i], "_tuning.pdf"), width = 10, height = 5)
}

