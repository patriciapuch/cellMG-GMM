# SIMULATION SETUP
library(mvtnorm)


###### GENERATE DATA #######
# covariance matrices
genData = function( p = 10,
                    n = 100,
                    N = 5,
                    cond = 100,
                    pi_diag = 0.75,
                    seed = 1,
                    m = 0,
                    cell_eps = 0.1,
                    gamma = 1,
                    row_eps = 0,
                    type = "ALYZCOR",
                    mu_type = "fr",  # or cs
                    csep = NULL) {


  # m: value for mu (eg 0 for all groups)
  # k: which eigenvector used for contamination (between 1/p - largest and 1-smallest)
  set.seed(seed)

  if(!type %in% c("ALYZCOR", "A0X")) stop("Wrong type for covariance structure!")

  # pi matrix
  # random:
  # pi_groups = matrix(abs(runif(N*N, min = 0, max = 1)), N, N)
  # not random:
  pi_groups = matrix(1, N, N)
  diag(pi_groups) = 0
  if(any(rowSums(pi_groups) == 0)) {
    pi_groups = diag(1, N)
  } else {
    pi_groups = (1-pi_diag)*pi_groups/(rowSums(pi_groups)) + diag(pi_diag, N)
  }

  # group and source assignments
  if(length(n) == 1) n = rep(n, N)
  if(length(n) != N) stop("n has not the size N!")
  groups = rep(1:N, times = n)
  source = rep(NA, sum(n))
  for(i in 1:N){
    start = ifelse(i == 1, 1, 1 + cumsum(n[1:(i-1)])[(i-1)])
    end = cumsum(n[1:i])[i]
    source[start:end] = sample(1:N, size = n[i], replace = T, prob = pi_groups[i, ])
  }

  Sigma = list()
  mu = list()

  for(i in 1:N){
    # construction of random matrices according to Agostinelli
    # without repetition of step 4 and 5, since we do not have group wise equivariance
    if (type == "ALYZCOR"){
      TT = matrix(rnorm(p*p), p, p)
      UU = eigen(TT %*% t(TT))$vectors
      dd = diag(sort(c(cond, 1, runif(p-2, min = 1, max = cond)), decreasing = T))
      S = UU %*% dd %*%  t(UU)

      # step 4
      R = cov2cor(S)

      # step 5
      U = eigen(R)$vectors
      d = eigen(R)$values
      d[1] = d[p]*cond
      R = U %*% diag(d) %*% t(U)

      while(sum(diag(R)) > 2*p | sum(diag(R)) < p*0.5){
        # step 4
        R = cov2cor(R)

        # step 5
        U = eigen(R)$vectors
        d = eigen(R)$values
        d[1] = d[p]*cond
        R = U %*% diag(d) %*% t(U)
      }


    }
    if(type == "A0X"){
      # correlation matrices but with random correlation between 0.9 and 0.5 and then make the same approach
      corr = runif(1, min = 0.5, max = 0.9)
      columns <- matrix(data = (seq_len(p)), nrow = p, ncol = p, byrow = TRUE)
      rows <- matrix(data = (seq_len(p)), nrow = p, ncol = p, byrow = FALSE)
      R <- matrix(-corr, nrow = p, ncol = p)
      R <- R^(abs(columns - rows))
    }
    Sigma = c(Sigma, list(R))

    if(mu_type == "fr"){  # if fixed repeated
      if (length(m) == 1 & is.numeric(m)){
        mu =  c(mu, list(rep(m, p)))
      } else {stop("For mu_type = 'fr' a numeric value need o be given.")}
    }
    if(mu_type == "cs"){
      if(is.null(csep)) stop("Please choose c vor c-separated clusters.")
      mu = c_sep_mean (Sigma, c = csep)
    }
  }

  # construct data
  Xclean = matrix(NA, nrow = sum(n), ncol = p)
  cont = matrix(1, sum(n), p)

  for(i in 1:N){
    # generate clean data per source
    ind = which(source == i)

    # clean data
    Xclean[ind, ] = mvtnorm::rmvnorm(n = length(ind), mean = mu[[i]], sigma = Sigma[[i]])
  }

  X_cont = Xclean

  for (i in 1:N){
    # generate contaminated data per group
    ind = which(groups == i)

    mhelp = matrix(FALSE, nrow = length(ind), ncol = p)
    ind_row = sample(1:dim(mhelp)[1], row_eps*dim(mhelp)[1])
    mhelp[ind_row, ] = TRUE

    mhelp = apply(X = mhelp,
                  MARGIN = 2,
                  FUN = function(x) {
                    y = x
                    if(length(ind_row) > 0 ){
                      ind_cells = sample((1:length(x))[-ind_row], size = cell_eps*length(x), replace = FALSE)
                    } else{
                      ind_cells = sample((1:length(x)), size = cell_eps*length(x), replace = FALSE)
                    }
                    y[ind_cells] = TRUE
                    y})
    contind = which(mhelp)

    cont[ind, ] = (!mhelp)*1

    for(l in ind){
      j = source[l]
      vars_cont =  which(cont[l,] == 0)

      if(length(vars_cont) > 0){
        # cellwise outliers
        eigsmall = eigen(Sigma[[j]][vars_cont, vars_cont, drop = FALSE])$vector
        eigsmall = eigsmall[, dim(eigsmall)[2]]
        md = c(t(eigsmall) %*%  solve(Sigma[[j]][vars_cont, vars_cont, drop = FALSE]) %*% c(eigsmall))
        X_cont[l, vars_cont] =  eigsmall * (gamma*sqrt(length(vars_cont))/sqrt(md)) + mu[[j]][vars_cont]
      }
    }
  }

  return(list(mu = mu,
              Sigma = Sigma,
              pi_groups = pi_groups,
              X = X_cont,
              groups = groups,
              source = source,
              W = cont
              ))
}

c_sep = function(mu1, mu2, Sigma1, Sigma2){

  p = length(mu1)

  a1 = sqrt(sum((mu1-mu2)^2))
  a2 = sqrt(p*eigen(Sigma1)$values[1] * eigen(Sigma2)$values[1])

  return(a1/a2)
}

c_sep_mean = function(Sigma, c = 1){

  N = length(Sigma)
  p = dim(Sigma[[1]])[1]

  mu = list(rep(0, p))
  for (k in 2:N){
    mu_mean = colMeans(do.call(rbind, mu))
    mu_random = rnorm(n = p)

    t = 0
    for(i in 1:k){
      f = function(t) c_sep(mu1 = mu[[i]], mu2 = mu_mean + t*(mu_random-mu_mean), Sigma1 = Sigma[[i]], Sigma2 = Sigma[[k]]) - c
      tmp <- try( uniroot(f = f, lower = 0, upper = 1000), silent = TRUE)
      if (!is(tmp, "try-error")) {
        t_c <- tmp$root
      } else {
        t_c = 0
      }
      if(t_c > t) t = t_c
    }
    mu = c(mu, list(mu_mean + t*(mu_random-mu_mean)))
  }
  mu

  return(mu)
}



###### OC METHOD #######
ollerer_croux = function(X){

  X = as.matrix(X)
  p = dim(X)[2]
  n = dim(X)[1]

  covs = matrix(NA, p, p)

  Q = apply(X, FUN = robustbase::Qn, MARGIN = 2)

  for(i in 1:p){
    for(j in 1:p){
      covs[i,j] = Q[i]*Q[j]*cor(X[, i], X[, j], method="kendall")
    }
  }
  return(covs)
}


##########################################################################################
# EVALUATION FUNCTIONS                                                                ####
##########################################################################################

eval_table = function(table){
  # positive is for outliers

  TP <- table[1, 1]  # True Positives
  FP <- table[2, 1]  # False Positives
  FN <- table[1, 2]  # False Negatives
  TN <- table[2, 2]  # True Negatives

  # Berechnung von Precision, Recall und F1-Score
  TPR <- TP / (TP + FN)  # True Positive Rate (Recall/Sensitivität)
  TNR <- TN / (TN + FP)  # True Negative Rate (Spezifität)
  FPR <- FP / (FP + TN)  # False Positive Rate
  FNR <- FN / (TP + FN)  # False Negative Rate

  precision <- TP / (TP + FP)   # percent of detected outliers is really outlying
  recall <- TP / (TP + FN)      # how many outliers are correctly identified
  f1 <- 2 * (precision * recall) / (precision + recall)

  return(c("TPR" = TPR,
           "TNR" = TNR,
           "FPR" = FPR,
           "FNR" = FNR,
           "precision" = precision,
           "recall" = recall,
           "f1" = f1))
}


best_ordering_covs = function(x, N){
  # Sort covariances and rest of output for mclust and cellGMM to get best performance
  # x: list element from simulation results including simulated data and output of method

  library(gtools)
  gruppen <- c(1:N)
  perms <- permutations(n = length(gruppen), r = length(gruppen), v = gruppen)
  perms = cbind(perms, NA)

  kl_min = Inf
  Sigma = x$out$Sigma
  mu = x$out$mu
  probs = x$out$probs
  for(i in 1:nrow(perms)){
    kl_cur = 0
    for(j in 1:N){
      kl_cur = kl_cur + KL(x$simulated$Sigma[[j]], x$out$Sigma[[perms[i, j]]])
    }

    # reorder if KL is lower as current solution
    if(kl_cur < kl_min){
      kl_min = kl_cur
      for(j in 1:N){
        Sigma[[j]] = x$out$Sigma[[perms[i, j]]]
        mu[[j]] = x$out$mu[[perms[i, j]]]
      }
      if(!any(is.na(x$out$probs)) & !any(is.na(x$out$pi_groups))){
        probs = x$out$probs[,perms[i,1:N]]
        pi_groups = x$out$pi_groups[perms[i,1:N]]
      }
    }
  }

  x$out$Sigma = Sigma
  x$out$mu = mu
  x$out$probs = probs
  x$out$pi_groups = NA

  return(x)
}


KL = function(S0, Shat){
  # Kullback-Leibler Divergence
  # S0: correct covariance matrix
  # Shat: estimated covariance matrix

  p = dim(S0)[1]
  kl = sum(diag(Shat %*% solve(S0))) - p - log(det(Shat %*% solve(S0)))
  return(kl)
}

MSE_mu = function(mu0, muhat){
  # MSE (mu)

  N = length(mu0)
  if(is.null(muhat)) {
    return(rep(NA, N))
  } else if(any(is.na(muhat))) {
    return(rep(NA, N))
  } else {
    return(sapply(1:N, function(x) mean((mu0[[x]] - muhat[[x]])^2)))
  }
}

MSE_pi = function(pi0, pihat){
  # MSE (pi)

  mse_pi = c(pi0-pihat)%*%c(pi0-pihat)/length(pi0)

  return(c(mse_pi))
}

MSE_probs= function(source, that, groups){
  # MSE (t_i)
  # source: real group/component vector with group number
  # that: estimated probability matrix
  # groups: pre-defined group with group number

  N = max(source)
  m = matrix(1:N, ncol = N, nrow = length(groups), byrow = TRUE)

  t_groups = as.numeric(groups == m)
  t_source = as.numeric(source == m)
  t_hat = c(that)

  index_noise = t_groups != t_source
  index_clean = t_groups == t_source


  MSE_all = c(t_source-t_hat)%*%c(t_source-t_hat)/length(t_hat)
  MSE_noise = c(t_source[index_noise]-t_hat[index_noise]) %*% c(t_source[index_noise]-t_hat[index_noise])/length(t_hat[index_noise])
  MSE_clean = c(t_source[index_clean]-t_hat[index_clean]) %*% c(t_source[index_clean]-t_hat[index_clean])/length(t_hat[index_clean])

  MSE_all_group = c(t_source-t_groups)%*%c(t_source-t_groups)/length(t_groups)
  MSE_noise_group = c(t_source[index_noise]-t_groups[index_noise]) %*% c(t_source[index_noise]-t_groups[index_noise])/length(t_groups[index_noise])
  MSE_clean_group = c(t_source[index_clean]-t_groups[index_clean]) %*% c(t_source[index_clean]-t_groups[index_clean])/length(t_groups[index_clean])

  return(c(MSE_all = MSE_all, MSE_noise = MSE_noise, MSE_clean = MSE_clean, MSE_all_group = MSE_all_group, MSE_noise_group = MSE_noise_group, MSE_clean_group = MSE_clean_group))
}

ARI = function(source, groups, class){
  # Adjusted Rand Index (from mclust)
  # source: real components
  # groups: pre-defined groups
  # class: estimated class by model

  ARI_group = mclust::adjustedRandIndex(as.numeric(source), as.numeric(groups))  # ARI pre-defined groups
  ARI_class = mclust::adjustedRandIndex(as.numeric(source), as.numeric(class))   # ARI estimated classes
  ARI_diff = ARI_class - ARI_group

  return(c(ARI_group = ARI_group, ARI_class = ARI_class, ARI_diff = ARI_diff))
}

class_out = function(Wreal, What, column = 0){
  # either for all or for certain variables
  if(column == 0)  column = 1:dim(Wreal)[2]
  wreal = factor(Wreal[, column], levels = c(0,1))
  what = factor(What[, column], levels = c(0,1))
  eval_table(table(wreal, what))
}


### HYPERPARAMETER TUNING ####
select_alpha_simulations = function(list_results, alpha_grid){

  obj = lapply(FUN = function (x) min(x$out$objvals, na.rm = T), X = list_results)
  obj = unlist(obj)
  ind = which.min(obj)
  #ind = which(diff(obj) >= 0)[1]
  if(is.na(ind)) ind = length(alpha_grid) # --> no decreasing step: ideal value alpha = 1

  alpha_select = alpha_grid[ind]

  # plot
  # plot(alpha_grid, obj, type = "l")
  # abline(v = alpha_select, col = "green")

  return(c(alpha = alpha_select, ind = ind, obj = obj))
}


select_alpha = function(list_results, alpha_grid){

  obj = lapply(FUN = function (x) min(x$objvals, na.rm = T), X = list_results)
  obj = unlist(obj)
  ind = which.min(obj)
  #ind = which(diff(obj) >= 0)[1]
  if(is.na(ind)) ind = length(alpha_grid) # --> no decreasing step: ideal value alpha = 1

  alpha_select = alpha_grid[ind]

  # plot
  # plot(alpha_grid, obj, type = "l")
  # abline(v = alpha_select, col = "green")

  return(c(alpha = alpha_select, ind = ind, obj = obj))
}
