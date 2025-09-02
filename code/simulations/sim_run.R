# RUN AND EVALUATE SIMULATIONS

# Load packages and code
library(foreach)
library(doParallel)
library(doSNOW)

# Load simulation setup code
path_sims ="sim_setup.R"
source(path_sims)

# Load cellGMM code originating from Github
source("cellGMM/cellGMM.R")
source("cellGMM/InitializationFunctions_cellGMM.R")
source("cellGMM/InternalFunctions_cellGMM.R")


# Construct parameter combinations
pnN = matrix(c(10, 100, 2,
               10, -2, 2,
               60, 40, 2,
               20, 30, 2,
               10, 100, 5),
             ncol = 3,
             byrow = TRUE)
cell_gamma = c(2, 6, 10)
pi_diag = c(0.75, 0.9)
csep = c(0, 0.5)
cell_eps = c(0.1)
seed = 1:100
method = c("cellGMM",
           "cellMGGMM",
           "ssMRCD",
           "MRCD",
           "Sample",
           "cellMCD",
           "ollerercroux",
           "mclust"
)
corr_type = c( "A0X", "ALYZCOR")
pars = expand.grid(cell_gamma = cell_gamma,
                   pi_diag = pi_diag,
                   csep = csep,
                   cell_eps = cell_eps,
                   seed = seed,
                   corr_type = corr_type,
                   method = method)
npars = dim(pars)[1]
pars = cbind(matrix(rep(pnN, each = npars),
                    byrow = F,
                    ncol = 3,
                    nrow = npars*dim(pnN)[1]),
             pars)
colnames(pars)[1:3] = c("p", "n", "N")
pars = pars[sort.int(pars$method, index.return = TRUE)$ix, ]
message(paste0("Number of Loops: ", dim(pars)[1]))


# Assign cores for workers
n_cores <- 100
message(paste0("Number of Cores: ",n_cores))
cl <- makeCluster(n_cores,  type = "SOCK")
registerDoSNOW(cl)

# Create the progress bar
pb <- txtProgressBar(min = 0, max = dim(pars)[1], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
progress_opts <- list(progress = progress)


# Run parallel loop
results <- foreach(i = 1:dim(pars)[1],
                   .packages = c("cellWise", "ssMRCD", "rrcov", "mclust", "tclust"),      # use ssMRCD Version 2.0.0
                   .options.snow = progress_opts,
                   .combine = c
                   ) %dopar% {

  a = Sys.time()
  out = list()
  out$pi_groups = NA
  out$probs = NA
  out$objvals = NA


  res = tryCatch({
    if(pars$n[i]>0) n = pars$n[i]

    # if n[i] = -100: unbalanced
    if(pars$n[i]==-2) n = c(100, 50)

    simulated = genData(p = pars$p[i],
                        n = n,
                        N = pars$N[i],
                        pi_diag = pars$pi_diag[i],
                        seed = pars$seed[i],
                        cell_eps = pars$cell_eps[i],
                        row_eps = 0,
                        gamma = pars$cell_gamma[i],
                        cond = 100,
                        type = pars$corr_type[i],
                        csep = pars$csep[i],
                        mu_type = "cs")

    start = Sys.time()
    if(pars$method[i] == "ollerercroux") {
      tmp = lapply(1:pars$N[i], function(x) ollerer_croux(X = simulated$X[simulated$groups == x, ]))
      out$Sigma = tmp
    }

    if(pars$method[i] == "mclust") {
      tmp = mclust::Mclust(data = simulated$X, G = pars$N[i])
      out$Sigma = lapply(1:pars$N[i], function(x) tmp$parameters$variance$sigma[, ,x])
      out$mu = lapply(1:pars$N[i], function(x) tmp$parameters$mean[, x])
    }

    if(pars$method[i] == "cellGMM") {
      # Hyperparameter setting:
      # From Supplement to Technometrics Paper for cellGMM:
      #   alpha_tclust =  2*alpha_true,
      #   alpha_1 = alpha_2 = alpha_true,
      #   alpha.A1 = alpha_true,
      #   alpha.A2 = 2*alpha.A1,
      #   nrep = 40,
      #   nstart = 10,
      #   niter = 10

      tmp = cellGMM(X = simulated$X,
                    G = pars$N[i],
                    tuning_param_init = data.frame(alpha_tclust  = 0.2,
                                                   alpha_1 = 0.1,
                                                   alpha_2 = 0.1,
                                                   alpha.A1 = 0.1,
                                                   alpha.A2 = 0.2,
                                                   nrep = 40,
                                                   nstart = 10,
                                                   niter = 10))
      out$Sigma = tmp$sigma
      out$mu = lapply(1:pars$N[i], function(x) tmp$mu[x, ])
      out$W = tmp$W
      out$pi_groups = tmp$pp
      out$probs = tmp$post
    }

    if(pars$method[i] == "cellMCD") {
      tmp = lapply(1:pars$N[i], function(x) cellWise::cellMCD(X = simulated$X[simulated$groups == x, ]))
      out$mu = lapply(tmp, function(x) x$mu)
      out$Sigma = lapply(tmp, function(x) x$S)
      out$W = do.call(rbind, lapply(tmp, function(x) x$W))
    }

    if(pars$method[i] == "MRCD") {
      tmp = lapply(1:pars$N[i], function(k) rrcov::CovMrcd(x = simulated$X[simulated$groups == k, ]))
      out$mu = lapply(tmp, function(x) x$center)
      out$Sigma = lapply(tmp, function(x) x$cov)
      out$hset = lapply(tmp, function(x) x$best)
    }

    if(pars$method[i] == "Sample") {
      out$mu = lapply(1:pars$N[i], function(k) colMeans(simulated$X[simulated$groups == k, ]))
      out$Sigma = lapply(1:pars$N[i], function(k) cov(simulated$X[simulated$groups == k, ]))
    }

    if(pars$method[i] == "ssMRCD") {
      if(pars$N[i] > 1){
        weights = matrix(1/(pars$N[i]-1), pars$N[i], pars$N[i])
        diag(weights) = 0
        tmp = ssMRCD::ssMRCD(X = simulated$X, groups = simulated$groups, lambda = 0.5, weights = weights)
        out$Sigma = tmp$Kcov
        out$mu = tmp$MRCDmu
        out$hset = tmp$hset
        out$rho = tmp$rho
        out$iter = tmp$numiter
      }
    }

    if(pars$method[i] == "cellMGGMM") {
      out = ssMRCD::cellMGGMM(X = simulated$X,
                              groups = simulated$groups,
                              alpha = 0.5,
                              hperc = 0.75,
                              nsteps = 500,
                              maxcond = 100)
    }

    end = Sys.time()
    res = list(out = out,
               simulated = simulated,
               time = difftime(end, start, units = c("secs")))
    }, error = function(cond){
      return(cond)
    })
  list(res)
  }


# Stop the parallel cluster after completion
close(pb)
stopCluster(cl)

# Save results
save(results, file = "simulations_result.RData")
save(pars, file = "simulations_parameters.RData")
##########################################################################################
# PLOTS                                                                               ####
##########################################################################################

# Load files
load(file = "simulations_result.RData")
load(file = "simulations_parameters.RData")

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


# Add scenario column
pars = mutate(pars,
              scenario = case_when(p == 10 & n == 100 & N == 2 ~ "balanced2",
                                   p == 10 & n == -2 ~ "unbalanced",
                                   p == 60 ~ "highdim",
                                   p == 20 ~ "mediumdim",
                                   p == 10 & n == 100 & N == 5 ~ "balanced5"))
head(pars)


# Check for non-successful runs
nonokay = which(sapply(results, function(x) "message" %in% names(x)))
okay = which(!1:dim(pars)[1] %in% nonokay)
if(length(nonokay) != 0) warning("Not all runs were successful!")

apply(pars[nonokay, ],
      MARGIN = 2,
      table)

table(sapply(nonokay, function(x) paste0(results[[x]])))                                  # only cellMCDlocal: Too many marginal outliers. or to many variables
apply(pars[1:nrow(pars) %in% nonokay & pars$method == "cellGMM", ],
      MARGIN = 2,
      table)                                                                              # and cellGMM: high-dimensional problems

# Check convergence of cellMG-GMM (positive increase of objective function)
tmp = sapply(which(pars$method == "cellMGGMM"),
             FUN = function(x) max(diff(results[[x]]$out$objvals)/min(results[[x]]$out$objvals, na.rm = T), na.rm = T))                                                                              # and cellGMM: high-dimensional problems
boxplot(tmp)  #-> computational inaccuracies

# Run time
time = unlist(sapply(results[okay], function(x) as.numeric(x$time)))
ggplot(pars[okay,] %>%
         cbind(time) %>%
         filter()) +
  geom_boxplot(aes(x = method,
                   y = time/60,
                   #col = scenario,
                   group = interaction(method, scenario))) +
  facet_grid(cols = vars(corr_type), rows = vars(scenario)) +
  scale_y_log10(breaks = c(0.0166666, 1, 60, 600), labels = c("1 s", "1 min", "1 h", "10 h")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.minor.y = element_blank()
        ) +
  labs(y ="", x = "")

# Use only successful runs
filtered_ind = okay

# Sort groups from cellGMM and mclust to get best performance
ind = which(pars$method %in% c("cellGMM", "mclust") & 1:nrow(pars) %in% okay)
for(i in ind){
  results[[i]] = best_ordering_covs(results[[i]], pars$N[i])
}


###### GENERATE PLOTS ######

variable_interest = matrix(c("sigma", "kl",
                             "w", "precision",
                             "w", "recall",
                             "w", "fscore",
                             "mu", "diff",
                             "probs", "mse"), ncol = 2, byrow = T)
pars_file = unique(pars[filtered_ind, c("p", "n", "N", "corr_type", "cell_eps", "scenario")]) %>%
  mutate(corr_type = paste(corr_type))

# For each scenario and covariance structure construct set of plots
for(i in 1:dim(pars_file)[1]){

  # Select indices of sub--scenarios
  indices_subset = which(do.call(paste, as.data.frame(pars[,c("p", "n", "N", "corr_type", "cell_eps")])) %in% do.call(paste, pars_file[i, 1:5]))
  indices_plot = indices_subset[indices_subset %in% filtered_ind]

  N_filter = pars_file$N[i]
  evals = matrix(NA, ncol = N_filter*3 + 19, nrow =length(indices_plot))
  colnames(evals) = c(paste0("sigma_kl", 1:N_filter),
                      "sigma_kl0",
                      paste0("mu_diff", 1:N_filter),
                      "mu_diff0",
                      paste0("probs_kl", 1:N_filter),
                      "probs_mse",
                      "probs_mae",
                      paste0("probs_", c("TPR","TNR","FPR", "FNR","precision", "recall","fscore")),  #7
                      paste0("w_", c("TPR","TNR","FPR", "FNR","precision", "recall","fscore")),  #7
                      "objval")  #1

  # Evaluate results for each scenario and covariance structure
  for(j in 1:length(indices_plot)){
    evals[j,] = unlist(eval_simul_simple(output = results[[indices_plot[j]]]$out,
                                         simulated = results[[indices_plot[j]]]$simulated))
  }

  # Clean evaluations
  eval_melted = evals %>%
    cbind(pars[indices_plot, ], .) %>%
    filter(cell_gamma  %in% c(2, 6, 10)) %>%
    melt(id.vars = colnames(pars)) %>%
    mutate(par = word(variable, 1, sep = "_"),
           metric = word(variable, 2, sep = "_")) %>%
    mutate(dist = as.numeric(gsub("[a-zA-Z]", "", metric)),
           metric = gsub("[0-9]", "", metric)) %>%
    filter(par %in% variable_interest[,1],
           metric %in% variable_interest[,2],
           (dist == 0 | is.na(dist)))

  # Plot evaluations
  for(pp in 1:dim(variable_interest)[1]){

    met = variable_interest[pp, 2]
    if(variable_interest[pp, 1] == "w") met = c("precision", "recall","fscore")

    dat_plot = eval_melted %>%
      filter(par == variable_interest[pp, 1], metric %in% met) %>%
      filter(!is.na(value)) %>%
      mutate(method = factor(case_when( method == "cellMGGMM" ~ "cellMG-GMM",
                                        method == "cellGMM" ~ "cellGMM",
                                        method == "ssMRCD" ~ "ssMRCD",
                                        method == "MRCD" ~ "MRCD",
                                        method == "Sample" ~ "sample",
                                        method == "cellMCD" ~ "cellMCD",
                                        method == "ollerercroux" ~ "OC",
                                        method == "mclust" ~ "mclust"),
                             levels =  c("cellMG-GMM", "cellMCD", "cellGMM", "OC", "MRCD", "ssMRCD", "mclust", "sample"),
                             ordered = TRUE))


    labels_fac = c("0" = "mu == 0",
                   "0.5" = "mu ~ ' varying'",
                   "0.9" = "pi[diag] == 0.9",
                   "0.75" = "pi[diag] == 0.75",
                   "recall" = "'Recall'",
                   "precision" = "'Precision'",
                   "fscore" = "'F1 Score'")

    colors = c("cellMG-GMM" = "#EE8866",
               "cellMCD" = "#EEDD88",
               "cellGMM" =  "#FFAABB",
               "OC" = "#BBCC33",
               "sample" = "#99DDFF",
               "ssMRCD" ="#44BB99",
               "MRCD" = "#AAAA00",
               "mclust" = "#77AADD")

    if(dim(dat_plot)[1] != 0 & variable_interest[pp, 1] == "sigma"){

       g = ggplot(data = dat_plot) +
        geom_boxplot(aes(x = as.factor(cell_gamma),
                         y = value,
                         fill = method,
                         group = interaction(method, cell_gamma)),
                     linewidth = 0.3,
                     fatten = 1.2,
                     outlier.size = 0.8) +
        facet_grid(rows = vars(csep),
                   cols = vars(pi_diag),
                   labeller = as_labeller(labels_fac, default = label_parsed),
                   scales = "free"
                   ) +
        labs(y = "KL-Divergence",
             x = expression(gamma[cell])
             ) +
        scale_y_log10() +
        scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
        scale_color_manual("", values = colors) +
        theme_bw(base_size = 16) +
        geom_vline(xintercept = seq(1.5, 2.5, 1), linewidth = 0.05, col = "grey") +
        theme(legend.position = "top",
              panel.grid = element_blank()
              ) +
         guides(fill = guide_legend(nrow = 1))

       g

       ggsave(g, file = paste0("simuls_", pars_file$scenario[i], "_", pars_file$corr_type[i], "_sigma.pdf"), width = 10, height = 5)
    }


    if(dim(dat_plot)[1] != 0 & variable_interest[pp, 1] == "mu"){

      g_mu = ggplot(data = dat_plot) +
        geom_boxplot(aes(x = as.factor(cell_gamma),
                         y = value,
                         fill = method,
                         group = interaction(method, cell_gamma)),
                     linewidth = 0.3,
                     fatten = 1.2,
                     outlier.size = 0.8) +
        facet_grid(rows = vars(csep),
                   cols = vars(pi_diag),
                   labeller = as_labeller(labels_fac, default = label_parsed),
                   scales = "free"
        ) +
        labs(y = expression(MSE(mu)),
             x = expression(gamma[cell])
             ) +
        scale_y_log10() +
        scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
        scale_color_manual("", values = colors) +
        theme_bw(base_size = 14) +
        geom_vline(xintercept = seq(1.5, 2.5, 1), linewidth = 0.05, col = "grey") +
        theme(legend.position = "top",
              panel.grid = element_blank()
        )

    }


    if(dim(dat_plot)[1] != 0 & variable_interest[pp, 1] == "probs"){

      g_probs = ggplot(data = dat_plot) +
        geom_boxplot(aes(x = as.factor(cell_gamma),
                         y = value,
                         fill = method,
                         group = interaction(method, cell_gamma)),
                     show.legend = FALSE,
                     linewidth = 0.3,
                     fatten = 1.2,
                     outlier.size = 0.8) +
        facet_grid(rows = vars(csep),
                   cols = vars(pi_diag),
                   labeller = as_labeller(labels_fac, default = label_parsed),
                   scales = "free"
        ) +
        labs(y = expression(MSE(pi)),
             x = expression(gamma[cell])
             ) +
        scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
        scale_color_manual("", values = colors) +
        theme_bw(base_size = 14) +
        theme(legend.position = "top",
              panel.grid = element_blank(),
              aspect.ratio = 3
        )
      g = grid.arrange(g_mu, g_probs, ncol=2, widths = c(2, 1))

      ggsave(g, file = paste0("simuls_", pars_file$scenario[i], "_", pars_file$corr_type[i], "_pimu.pdf"), width = 10, height = 6)
    }


    if(dim(dat_plot)[1] != 0 & variable_interest[pp, 1] == "w"){

      g = ggplot(data = dat_plot) +
        geom_boxplot(aes(x = as.factor(cell_gamma),
                         y = value,
                         fill = method,
                         group = interaction(method, cell_gamma)),
                     linewidth = 0.3,
                     fatten = 1.2,
                     outlier.size = 0.8) +
        facet_nested(csep ~ metric + pi_diag,
                     labeller = as_labeller(labels_fac, default = label_parsed)) +
        labs(y = "",
             x = expression(gamma[cell])
             ) +
        ylim(c(0,1)) +
        scale_fill_manual("", values = sapply(colors, function(x) scales::alpha(x, 0.75))) +
        scale_color_manual("", values = colors) +
        theme_bw(base_size = 16) +
        geom_vline(xintercept = seq(1.5, 2.5, 1), linewidth = 0.05, col = "grey") +
        theme(legend.position = "top",
              panel.grid = element_blank()
        ) +
        guides(fill = guide_legend(nrow = 1))
      g


      ggsave(g, file = paste0("simuls_", pars_file$scenario[i], "_", pars_file$corr_type[i], "_W.pdf"), width = 10, height = 5)
    }
  }

  # Create table for successful runs
  number_reps = eval_melted %>%
    select(cell_gamma, pi_diag, csep, seed, method) %>%
    unique() %>%
    group_by (cell_gamma, pi_diag, csep, method) %>%
    rename(mu = csep) %>%
    mutate(mu = ifelse(mu == 0, "0", "varying"),
           cell_gamma = as.character(cell_gamma)) %>%
    summarise(n = n()) %>%
    data.frame

  print.xtable(xtable(number_reps),
               include.rownames = FALSE,
               file = paste0("simuls_", pars_file$scenario[i], "_", pars_file$corr_type[i], "_nreps.txt"))
}
