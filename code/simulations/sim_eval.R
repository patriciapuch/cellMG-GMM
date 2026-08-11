##########################################################################################
# PLOTS                                                                               ####
##########################################################################################

#### SETUP ####

# Load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(gridExtra)
library(knitr)
library(kableExtra)
library(patchwork)
library(ggh4x)
library(xtable)
library(cowplot)

# Load code
source("sim_setup.R")

# Load files
load(file = "sim_results_fullALYZCOR.RData")
load(file = "sim_parameters_fullALYZCOR.RData")

# Clean and expand data frame
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
              method_factor = factor(case_when(  method == "cellMGGMM" ~ "cellMG-GMM",
                                                 method == "cellGMM" ~ "cellGMM",
                                                 method == "MGGMM" ~ "MG-GMM",
                                                 method == "ssMRCD" ~ "ssMRCD",
                                                 method == "MRCD" ~ "MRCD",
                                                 method == "Sample" ~ "sample",
                                                 method == "cellMCD" ~ "cellMCD",
                                                 method == "ollerercroux" ~ "OC",
                                                 method == "mclust" ~ "mclust",
                                                 method == "rmda" ~"rmda"),
                                      levels =  c("cellMG-GMM", "MG-GMM", "cellGMM", "cellMCD", "OC", "ssMRCD", "MRCD", "rmda", "mclust", "sample"),
                                      ordered = TRUE),
              contaminated = factor(case_when(cell_gamma == 0 ~ "No contamination",
                                       cell_gamma == 6 ~ "Contamination",
                                       TRUE ~ NA_character_), levels = c("No contamination", "Contamination"), ordered = TRUE) )


# Set colour/label scheme
colors = c("cellMG-GMM" = "#AA3377",
           "MG-GMM" = "#FFAABB",
           "cellGMM" =  "#EE8866",
           "cellMCD" = "#EEDD88",
           "OC" = "#BBCC33",
           "ssMRCD" ="#44BB99",
           "MRCD" = "#AAAA00",
           "rmda" = "#99DDFF",
           "mclust" = "#77AADD",
           "sample" = "#332288"
           )
labels_facets = c( "0" = "mu == 0",
                   "0.5" = "mu ~ ' varying'",
                   "0.9" = "pi[diag] == 0.9",
                   "1" = "pi[diag] == 1",
                   "0.75" = "pi[diag] == 0.75",
                   "recall" = "'Recall'",
                   "precision" = "'Precision'",
                   "fscore" = "'F1 Score'",
                   "MSE_all" = "'All observation'",
                   "MSE_noise" = "'Mislabelled observation'",
                   "MSE_clean" = "'Correctly labelled observation'",
                   "ARI_class" = "'ARI'",
                   "ARI_diff" = "'ARI (difference)'",
                   "ARI_group" = "'ARI (pre-defined)'",
                   "No contamination" = "'No contamination'",
                   "Contamination" = "'Contamination'")




#### TEST FOR FAILED RUNS ####
# Check for non-successful runs, dependence on parameters and error messages
nonokay = which(sapply(results, function(x) "message" %in% names(x) | is.null(x$out)))
okay = which(!1:dim(pars)[1] %in% nonokay)
if(length(nonokay) != 0) warning("Not all runs were successful!")

apply(pars[nonokay, ],
      MARGIN = 2,
      table)

apply(pars[1:nrow(pars) %in% nonokay & pars$method == "cellGMM", ],
      MARGIN = 2,
      table)

apply(pars[1:nrow(pars) %in% nonokay & pars$method == "cellMCD", ],
      MARGIN = 2,
      table)

apply(pars[1:nrow(pars) %in% nonokay & pars$method == "rmda", ],
      MARGIN = 2,
      table)

nonokay_message = which(sapply(results, function(x) "message" %in% names(x)))
table(sapply(which(1:nrow(pars) %in% nonokay_message & pars$method == "rmda"), function(x) paste0(results[[x]])))


# Check convergence of cellMG-GMM (positive increase of objective function)
tmp = sapply(which(pars$method == "cellMGGMM"),
             FUN = function(x) max(diff(results[[x]]$out$objvals)/min(results[[x]]$out$objvals, na.rm = T), na.rm = T))
boxplot(tmp)  #-> computational inaccuracies
tmp = sapply(which(pars$method == "cellMGGMM"),
             FUN = function(x) max(diff(results[[x]]$out$objvals), na.rm = T))
boxplot(tmp)  #-> computational inaccuracies



# filter for sucessful runs
results_success = results[okay]
pars_success = pars[okay, ]




#### SORTING FOR CLUSTERING METHODS ####
# Sort groups from cellGMM and mclust to get best performance
ind = which(pars_success$method %in% c("cellGMM", "mclust"))
for(i in ind){
  results_success[[i]] = best_ordering_covs(results_success[[i]], pars_success$N[i])
}




#### RUNTIME PLOT ####
time = unlist(sapply(results_success, function(x) as.numeric(x$time)))
pars_success = pars_success %>% mutate(time = time)

gg_time = ggplot(pars_success %>% filter(cell_gamma == 6,
                                 pi_diag == 0.75,
                                 csep == 0.5)) +
  geom_boxplot(aes(x = method_factor,
                   y = time/60,
                   fill = method_factor,
                   group = interaction(method_factor, scenario)),
               linewidth = 0.3,
               outlier.size = 0.4) +
  facet_wrap(vars(scenario_nr), ncol = 5) +
  scale_y_log10(breaks = c(0.0166666, 1, 60, 600), labels = c("1 s", "1 min", "1 h", "10 h")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none"
  ) +
  labs(y ="", x = "") +
  scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.85)))
gg_time

ggsave(gg_time,
  file = "simulation_runtime.pdf",
  width = 8, height = 3.5
)

#### Make plots for each scenario and covariance structure ####

combinations = unique(pars_success[, c("corr_type", "scenario")])
for(i in 1:nrow(combinations)){

  index_subset = which(pars_success$corr_type == combinations[i, 1] & pars_success$scenario == combinations[i, 2] &
                         pars_success$cell_gamma %in% c(0,6) & pars_success$method != "Sample")
  N_filter = 2
  if(combinations[i, 2] == "balanced5") N_filter = 5


  #### Covariance KL ####
  kl_plot = rep(NA, length(index_subset))
  for(j in 1:length(index_subset)){
    jj = index_subset[j]
      kl_sigma = sapply(1:length(results_success[[jj]]$simulated$Sigma), FUN = function(x) KL(results_success[[jj]]$simulated$Sigma[[x]], results_success[[jj]]$out$Sigma[[x]]))
      kl_plot[j] = mean(kl_sigma)

  }
  dat_plot = mutate(pars_success[index_subset, ], kl = kl_plot)

  g_KL = ggplot(data = dat_plot) +
    geom_boxplot(aes(x = factor(pi_diag, levels = c(1,0.9, 0.75)),
                     y = kl,
                     fill = method_factor,
                     group = interaction(cell_gamma, method_factor, pi_diag)),
                 linewidth = 0.3, median.linewidth = 0.3,
                 outlier.size = 0.1) +
    facet_grid(rows = vars(csep),
               cols = vars(contaminated),
               labeller = as_labeller(labels_facets,
                                      default = label_parsed),
               scales = "free"
    ) +
    labs(y = "KL-Divergence",
         x = expression(pi[diag]) # expression(gamma[cell])
    ) +
    scale_y_log10() +
    scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.85))) +
    theme_bw(base_size = 15) +
    geom_vline(xintercept = seq(1.5, 2.5, 1),
               linewidth = 0.05,
               col = "grey") +
    theme(panel.grid = element_blank(),
          #legend.position = "top",
          plot.margin = margin(1, 1, 1, 1)
    ) +
    guides(fill = guide_legend(ncol = 1))
  g_KL

  ggsave(g_KL, file = paste0("simuls_", combinations$scenario[i], "_", combinations$corr_type[i], "_sigma.pdf"), width = 10, height = 4.5)



  #### Mean and pi ####
  mse_mu_plot = rep(NA, length(index_subset))
  mse_pi_plot = rep(NA, length(index_subset))

  for(j in 1:length(index_subset)){
    jj = index_subset[j]
    mse_mu = sapply(1:length(results_success[[jj]]$simulated$Sigma),
                    FUN = function(x) MSE_mu(results_success[[jj]]$simulated$mu[[x]], results_success[[jj]]$out$mu[[x]]))
    mse_mu_plot[j] = mean(mse_mu)

    mse_pi_plot[j] = if(any(!is.na(results_success[[jj]]$out$pi_groups))){
                        MSE_pi(pi0 = results_success[[jj]]$simulated$pi_groups,
                               pihat = results_success[[jj]]$out$pi_groups)
                      } else {NA}
  }
  dat_plot = mutate(pars_success[index_subset, ],
                    mse_mu = mse_mu_plot,
                    mse_pi = mse_pi_plot)

  g_mu = ggplot(data = dat_plot) +
    geom_boxplot(aes(x = factor(pi_diag, levels = c(1,0.9, 0.75)), # as.factor(cell_gamma),
                     y = mse_mu,
                     fill = method_factor,
                     group = interaction(method_factor, cell_gamma, pi_diag)),
                 linewidth = 0.3, median.linewidth = 0.3,
                 outlier.size = 0.1) +
    facet_grid(rows = vars(csep),
               cols = vars(contaminated), #cols = vars(pi_diag),
               labeller = as_labeller(labels_facets, default = label_parsed),
               scales = "free"
    ) +
    labs(y = expression(MSE(mu)),
         x = expression(pi[diag])#expression(gamma[cell])
    ) +
    scale_y_log10() +
    scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
    theme_bw(base_size = 13) +
    geom_vline(xintercept = seq(1.5, 2.5, 1), linewidth = 0.05, col = "grey") +
    theme(#legend.position = "top",
          panel.grid = element_blank(),
          plot.margin = margin(1, 1, 1, 1)
    )
  g_mu


  g_pi = ggplot(data = dat_plot) +
    geom_boxplot(aes(x = factor(pi_diag, levels = c(1,0.9, 0.75)), #as.factor(cell_gamma),
                     y = mse_pi,
                     fill = method_factor,
                     group = interaction(method_factor, cell_gamma,pi_diag)),
                 #show.legend = FALSE,
                 linewidth = 0.3, median.linewidth = 0.3,
                 outlier.size = 0.1) +
    facet_grid(rows = vars(csep),
               cols = vars(contaminated),# cols = vars(pi_diag),
               labeller = as_labeller(labels_facets, default = label_parsed),
               scales = "free"
    ) +
    labs(y = expression(MSE(pi)),
         x = expression(pi[diag])# expression(gamma[cell])
    ) +
    scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top",
          panel.grid = element_blank(),
          plot.margin = margin(1, 1, 1, 20)
    ) +
    geom_vline(xintercept = seq(1.5, 2.5, 1), linewidth = 0.05, col = "grey") +
    guides(fill = guide_legend(ncol = 1))
  g_pi


  g = grid.arrange(g_mu, g_pi, ncol=2, widths = c(2, 1))
  g
  ggsave(g, file = paste0("simuls_", combinations$scenario[i], "_", combinations$corr_type[i], "_pimu.pdf"), width = 10, height = 5)



  #### W ####
  w_eval = matrix(NA,
                  ncol = 7,
                  nrow = length(index_subset),
                  dimnames = list(NULL, c("TPR","TNR","FPR", "FNR","precision", "recall","fscore")))
  for(j in 1:length(index_subset)){
    jj = index_subset[j]
    w_eval[j, ] = if(!is.null(results_success[[jj]]$out$W)){
      class_out(Wreal = results_success[[jj]]$simulated$W,
                What = results_success[[jj]]$out$W)
    } else {rep(NA, 7)}
  }
  dat_plot = cbind(pars_success[index_subset, ], w_eval) %>%
    tidyr::pivot_longer(cols = TPR:fscore) %>%
    filter(name %in% c("precision", "recall","fscore"),
           method_factor != "MG-GMM",
           cell_gamma != 0)

  g_W = ggplot(data = dat_plot) +
    geom_boxplot(aes(x = factor(pi_diag, levels = c(1,0.9, 0.75)),# as.factor(cell_gamma),
                     y = value,
                     fill = method_factor,
                     group = interaction(method_factor, cell_gamma, pi_diag)),
                 linewidth = 0.3, median.linewidth = 0.3,
                 outlier.size = 0.1) +
    facet_nested(csep ~ name,
                 labeller = as_labeller(labels_facets, default = label_parsed)) +
    labs(y = "",
         x = expression(pi[diag])#expression(gamma[cell])
    ) +
    ylim(c(0,1)) +
    scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
    theme_bw(base_size = 14) +
    geom_vline(xintercept = c(1.5, 2.5),
               linewidth = 0.05,
               col = "grey") +
    theme(legend.position = "top",
          panel.grid = element_blank(),
          plot.margin = margin(1, 1, 1, 1)
    ) +
    guides(fill = guide_legend(nrow = 1))
  g_W

  ggsave(g_W, file = paste0("simuls_", combinations$scenario[i], "_", combinations$corr_type[i], "_W.pdf"), width = 10, height = 4)




  #### Class Probabilities ####

  t_eval = matrix(NA,
                  ncol = 6,
                  nrow = length(index_subset),
                  dimnames = list(NULL, c("ARI_group", "ARI_class", "ARI_diff", "MSE_all", "MSE_noise", "MSE_clean")))
  for(j in 1:length(index_subset)){
    jj = index_subset[j]
    if(!is.null(results_success[[jj]]$out$probs) & any(!is.na(results_success[[jj]]$out$probs))){
      class = apply(X = results_success[[jj]]$out$probs, FUN = which.max, MARGIN = 1)
      t_eval[j, 1:3] = ARI(source = results_success[[jj]]$simulated$source,
                           groups = results_success[[jj]]$simulated$groups,
                           class = class)
      t_eval[j, 4:6] = MSE_probs(source = results_success[[jj]]$simulated$source,
                                 that = results_success[[jj]]$out$probs,
                                 groups = results_success[[jj]]$simulated$groups)[1:3]
    } else {
      t_eval[j, ] = rep(NA, 6)
    }
  }

  dat_plot = cbind(pars_success[index_subset, ], t_eval) %>%
    tidyr::pivot_longer(cols = ARI_group:MSE_clean) %>%
    filter(name %in% c("ARI_group", "ARI_class", "ARI_diff", "MSE_all", "MSE_noise", "MSE_clean"))

  g_that_ari = ggplot(data = dat_plot %>%
                    filter(name %in% c( "ARI_class"))) +
    geom_boxplot(aes(x = factor(pi_diag, levels = c(1,0.9, 0.75)),#as.factor(cell_gamma),
                     y = value,
                     fill = method_factor,
                     group = interaction(method_factor, cell_gamma, pi_diag)),
                 linewidth = 0.3, median.linewidth = 0.3,
                 outlier.size = 0.1) +
    facet_nested(csep ~  contaminated,
                 labeller = as_labeller(labels_facets,
                                        default = label_parsed)) +
    labs(y = "ARI",
         x = expression(pi[diag])#expression(gamma[cell])
    ) +
    ylim(c(-1,1)) +
    scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
    theme_bw(base_size = 14) +
    geom_vline(xintercept = seq(1.5, 2.5, 1), linewidth = 0.05, col = "grey") +
    theme(legend.position = "top",
          panel.grid = element_blank(),
          plot.margin = margin(1, 1, 1, 1)
    ) +
    guides(fill = guide_legend(nrow = 1))
  g_that_ari
  ggsave(g_that_ari, file = paste0("simuls_", combinations$scenario[i], "_", combinations$corr_type[i], "_that_ari.pdf"), width = 10, height = 5)


  g_that_mse = ggplot(data = dat_plot %>%
                        filter(name %in% c("MSE_all", "MSE_noise"))) +
    geom_boxplot(aes(x = factor(pi_diag, levels = c(1,0.9, 0.75)),#as.factor(cell_gamma),
                     y = value,
                     fill = method_factor,
                     group = interaction(method_factor, cell_gamma,pi_diag)),
                 linewidth = 0.3, median.linewidth = 0.3,
                 outlier.size = 0.1) +
    # geom_hline(aes(yintercept = 0.25)) +
    # geom_hline(aes(yintercept = 0.1)) +
    facet_nested(csep ~ contaminated + name ,
                 labeller = as_labeller(labels_facets,
                                        default = label_parsed)) +
    labs(y = "MSE(t)",
         x = expression(pi[diag])#expression(gamma[cell])
    ) +
    scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
    theme_bw(base_size = 14) +
    geom_vline(xintercept = seq(1.5, 2.5, 1), linewidth = 0.05, col = "grey") +
    theme(legend.position = "top",
          panel.grid = element_blank(),
          plot.margin = margin(1, 1, 1, 1)
    ) +
    guides(fill = guide_legend(nrow = 1))
  g_that_mse

  ggsave(g_that_mse, file = paste0("simuls_", combinations$scenario[i], "_", combinations$corr_type[i], "_that_mse.pdf"), width = 10, height = 5)

}


# Create table for successful runs
number_reps = pars_success %>%
  select(cell_gamma, pi_diag, csep, seed, scenario, method) %>%
  filter(cell_gamma  != 10) %>%
  unique() %>%
  group_by (cell_gamma, pi_diag, csep, scenario, method) %>%
  rename(mu = csep) %>%
  mutate(mu = ifelse(mu == 0, "0", "varying"),
         cell_gamma = as.character(cell_gamma)) %>%
  summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = c(method), values_from = n, values_fill = 0) %>%
  select(where(~mean(.x == 100) != 1)) %>%
  arrange(scenario, as.numeric(cell_gamma)) %>%
  data.frame
number_reps

print.xtable(xtable(number_reps),
             include.rownames = FALSE,
             file = paste0("simuls_nreps_all_",  combinations$corr_type[i], ".txt"))



