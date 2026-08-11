# RUN AND EVALUATE SIMULATIONS

# Load packages and code
library(foreach)
library(doParallel)
library(doSNOW)

# Load simulation setup code
source("sim_setup.R")

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
cell_gamma = c(0, 2, 6, 10)  #0: no contamination
pi_diag = c(0.75, 0.9, 1)
csep = c(0, 0.5)
seed = 1:100
method = c("cellGMM",
           "cellMGGMM",
           "MGGMM",
           "ssMRCD",
           "MRCD",
           "Sample",
           "cellMCD",
           "ollerercroux",
           "mclust",
           "rmda"
)
#corr_type = c("A0X")
corr_type = c("ALYZCOR")
pars = expand.grid(cell_gamma = cell_gamma,
                   pi_diag = pi_diag,
                   csep = csep,
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
n_cores <- 4
message(paste0("Number of Cores: ",n_cores))
cl <- makeCluster(n_cores,  type = "SOCK")
registerDoSNOW(cl)

# Create the progress bar
pb <- txtProgressBar(min = 0, max = dim(pars)[1], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
progress_opts <- list(progress = progress)


# Run parallel loop
results <- foreach(i = 1:dim(pars)[1],
                   .packages = c("cellWise", "ssMRCD", "rrcov", "mclust", "tclust", "robustDA"),      # use ssMRCD Version 2.0.0
                   .options.snow = progress_opts,
                   .combine = c
                   ) %dopar% {

  out = list()
  out$pi_groups = NA
  out$probs = NA
  out$objvals = NA

  # set contamination
  if(pars$cell_gamma[i] == 0) {
    cell_eps = 0
  } else {
    cell_eps = 0.1
  }


  res = tryCatch({
    if(pars$n[i]>0) n = pars$n[i]

    if(pars$n[i]==-2) n = c(100, 50)

    simulated = genData(p = pars$p[i],
                        n = n,
                        N = pars$N[i],
                        pi_diag = pars$pi_diag[i],
                        seed = pars$seed[i],
                        cell_eps = cell_eps,
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
      out$probs = tmp$z
    }

    if(pars$method[i] == "rmda") {
      tmp = robustDA::rmda(X = simulated$X, cls = simulated$groups, K = pars$N[i], model = "VVV")
      out$Sigma = lapply(1:pars$N[i], function(x) tmp$prms$parameters$variance$sigma[, ,x])
      out$mu = lapply(1:pars$N[i], function(x) tmp$prms$parameters$mean[, x])
      out$pi_groups = tmp$R
      out$probs = tmp$prms$z
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
                              alpha = 0.75,
                              hperc = 0.75,
                              nsteps = 500,
                              maxcond = 100)
    }

    if(pars$method[i] == "MGGMM") {
      out = ssMRCD::cellMGGMM(X = simulated$X,
                              groups = simulated$groups,
                              alpha = 0.75,
                              hperc = 1,
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
save(results, file = paste0("sim_result", corr_type, ".RData"))
save(pars, file = paste0("sim_parameters", corr_type, ".RData"))
