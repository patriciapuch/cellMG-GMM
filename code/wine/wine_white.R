# WHITE WINE DATA (UCI Machine Learning Repository)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ssMRCD)
library(ellipse)
library(geomtextpath)

source("../simulations/sim_setup.R")
source("../simulations/cellGMM/cellGMM.R")
source("../simulations/cellGMM/InternalFunctions_cellGMM.R")
source("../simulations/cellGMM/InitializationFunctions_cellGMM.R")

# read in data
wine_white = read_delim("winequality-white.csv",
                        delim = ";",
                        escape_double = FALSE,
                        trim_ws = TRUE)

# construct groups
table(wine_white$quality)
groups_white = wine_white$quality
groups_white[groups_white == 9] = 7
groups_white[groups_white == 8] = 7
groups_white[groups_white == 3] = 5
groups_white[groups_white == 4] = 5
groups_white = groups_white - 4
table(groups_white)

# remove quality column
wine_white = wine_white[, -12]

# use robust transformation to normality
wine_white = cellWise::transfo(wine_white)$Y

##########################################################################################
# Hyper parameter tuning for alpha (time: ~5 min)
alpha_grid = seq(0.5, 1, 0.05)
l = vector(mode = "list", length = length(alpha_grid))
for(i in 1:length(alpha_grid)){
  l[[i]] = ssMRCD::cellMGGMM(X = wine_white,
                             groups = groups_white,
                             nsteps = 100,
                             alpha = alpha_grid[i],
                             maxcond = 100)
}
save(l, file ="white_wine_parametertuning_tmp.RData")g

out_tune = select_alpha(list_results = l, alpha_grid = alpha_grid)
plot(alpha_grid, out_tune[-c(1,2)], type = "l") # alpha: 0.65
##########################################################################################

# optimal model
out = l[[out_tune["ind"]]]

out$rho
# mixture probabilities
cat("Pi (in %):\n")
round(out$pi_groups*100, 2)

# percentage of outliers
cat("% Outliers per group and variable:\n")
round(sapply(1:length(unique(groups_white)), function(x) colMeans(1-out$W[groups_white == x, ]))*100, 2)

# calculate residuals
res = ssMRCD::residuals_mggmm(X = wine_white,
                     groups = groups_white,
                     Sigma = out$Sigma,
                     mu = out$mu,
                     probs = out$probs,
                     W = out$W,
                     set_to_zero = TRUE)

outs_wine = wine_white
outs_wine[res == 0] = NA

# parallel coordinate plot - preparation
dat_all = rbind(cbind(data.frame(wine_white),
                      kind = "data",
                      groups = groups_white,
                      assign = as.numeric(apply(out$probs, 1, which.max))),
                cbind(data.frame(outs_wine),
                      kind = "outs",
                      groups = groups_white,
                      assign = as.numeric(apply(out$probs, 1, which.max))),
                cbind(data.frame(rbind(out$mu[[1]],
                                       out$mu[[2]],
                                       out$mu[[3]])),
                      kind = "mu",
                      groups = rep(1:3, times = 3),
                      assign = rep(1:3, each = 3)),
                cbind(data.frame(sqrt(rbind(diag(out$Sigma[[1]]),
                                            diag(out$Sigma[[2]]),
                                            diag(out$Sigma[[3]])))),
                      kind = "S",
                      groups = rep(1:3, times = 3),
                      assign = rep(1:3, each = 3)),
                cbind(data.frame(rbind(out$mu[[1]] + sqrt(diag(out$Sigma[[1]])),
                                       out$mu[[2]] + sqrt(diag(out$Sigma[[2]])),
                                       out$mu[[3]] + sqrt(diag(out$Sigma[[3]])))),
                      kind = "upper",
                      groups = rep(1:3, times = 3),
                      assign = rep(1:3, each = 3)),
                cbind(data.frame(rbind(out$mu[[1]] - sqrt(diag(out$Sigma[[1]])),
                                       out$mu[[2]] - sqrt(diag(out$Sigma[[2]])),
                                       out$mu[[3]] - sqrt(diag(out$Sigma[[3]])))),
                      kind = "lower",
                      groups = rep(1:3, times = 3),
                      assign = rep(1:3, each = 3)))

dat_all$id = 1:dim(dat_all)[1]


ord = c("fixed acidity","volatile acidity","citric acid","residual sugar","chlorides",
        "free sulfur dioxide","total sulfur dioxide", "density","pH","sulphates","alcohol")
datplot = tidyr::pivot_longer(dat_all, "fixed.acidity":"alcohol") %>%
  filter(kind != "S") %>%
  mutate(name = gsub("\\.", " ", name)) %>%
  mutate( name = factor(name, ordered = TRUE, level = ord),
         groups = factor(case_when(groups == 1 ~ "Low Quality (Experts)",
                            groups == 2 ~ "Middle Quality (Experts)",
                            groups == 3 ~ "High Quality (Experts)",
                            TRUE ~ NA_character_),
                         ordered = TRUE,
                         level = c("Low Quality (Experts)", "Middle Quality (Experts)", "High Quality (Experts)")),
         assign = factor(case_when(assign == 1 ~ "Low Quality \n(Model)",
                            assign == 2 ~ "Middle Quality \n(Model)",
                            assign == 3 ~ "High Quality \n(Model)",
                            TRUE ~ NA_character_),
                         ordered = TRUE,
                         level = c("Low Quality \n(Model)", "Middle Quality \n(Model)", "High Quality \n(Model)")))

# parallel coordinate plot - plot
g_par = ggplot() +
  geom_line(data = datplot %>% filter(kind == "data"),
            aes(x = name,
                y = value,
                col = as.factor(groups),
                group = id,
                alpha = interaction(assign, groups)
                ),
            linewidth = 0.5) +
  geom_line(data = datplot %>% filter(kind == "mu"),
            aes(x = name,
                y = value,
                group = id),
            col = "black",
            linewidth = 0.75) +
  geom_errorbar(data = datplot %>%
                filter(kind != "data", kind != "outs") %>%
                select(-id) %>%
                  tidyr::pivot_wider(names_from = kind, values_from = value),
              aes(ymin = lower,
                  ymax = upper,
                  x = name
                  ),
              linewidth = 0.75,
              width = 0.2
              ) +
  theme_bw(base_size = 16) +
  scale_size_manual(values = c(0.4, 1,1, 1))+
  scale_color_manual(values = c("darkblue", "darkgreen", "darkred")) +
  scale_alpha_manual(values = c(0.02, 0.05, 0.05,
                                0.05, 0.02, 0.03,
                                0.1, 0.04, 0.02)) +
  scale_y_continuous(transform = "pseudo_log",
                     breaks = c(0, 10, 30)) +
  facet_grid(cols = vars(groups),
             rows = vars(assign)) +
  labs(x = "", y = "", title = "") +
  theme(legend.position = "top",
        title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(0, 0, 0, 20)) +
  guides(alpha = "none",
         col = "none")

g_par
ggsave(g_par, file = "wine_white_pcp_065.pdf", width = 10, height = 7)


####################################
## Ellipse Plot

# duration: 24h
out_cellGMM = cellGMM(X = wine_white,
                      G = 3,
                      tuning_param_init = data.frame(alpha_tclust  = 0.2,
                                                     alpha_1 = 0.1,
                                                     alpha_2 = 0.1,
                                                     alpha.A1 = 0.1,
                                                     alpha.A2 = 0.2,
                                                     nrep = 10,
                                                     nstart = 5,
                                                     niter = 10))
save(out_cellGMM, file = "wine_white_cellGMM.RData")


var1 = "alcohol"
var2 = "residual.sugar"
which.var = c(which(colnames(wine_white) == var1), which(colnames(wine_white) == var2))
data_plot = wine_white

ellipse_data = tibble(!!var1 := numeric(),
                      !!var2 := numeric(),
                      method = character(),
                      group = numeric())
for(i in 1:3){

  # cellMG-GMM
  cov = out$Sigma[[i]][which.var, which.var]
  mean = out$mu[[i]][which.var]
  ellipse_data = rbind(ellipse_data,
                       as.data.frame(ellipse(cov, centre = mean, level = 0.95)) %>%
                         mutate(method = "cellMG-GMM", group = i))

  # cellMG-GMM (alpha = 1)
  cov = l[[length(l)]]$Sigma[[i]][which.var, which.var]
  mean = l[[length(l)]]$mu[[i]][which.var]
  ellipse_data <- rbind(ellipse_data,
                        as.data.frame(ellipse(cov, centre = mean, level = 0.95)) %>%
                          mutate(method = "cellMG-GMM (alpha = 1)", group = i))

  # cellGMM
  reorderind = c(3, 1, 2)
  cov = out_cellGMM$sigma[[i]][which.var, which.var]
  mean = out_cellGMM$mu[i,][which.var]
  ellipse_data <- rbind(ellipse_data,
                        as.data.frame(ellipse(cov, centre = mean, level = 0.95)) %>%
                          rename(!!var1 := x, !!var2 := y) %>%
                          mutate(method = "cellGMM", group = reorderind[i]))

}

data_plot = cbind(data.frame(rbind(wine_white,wine_white, wine_white)),
      group = c(groups_white, groups_white, groups_white),
      class= c(groups_white, as.numeric(out$class), reorderind[out_cellGMM$label]),
      method = rep(c("cellMG-GMM (alpha = 1)", "cellMG-GMM", "cellGMM"), each = length(groups_white)))

cols = c("darkblue", "darkgreen", "darkred")
quality_labels <- c(
  "1"                      = "Low~Quality~(Experts)",
  "2"                      = "Middle~Quality~(Experts)",
  "3"                      = "High~Quality~(Experts)",
  "cellMG-GMM (alpha = 1)" = 'atop("cellMG-GMM", (alpha == 1))',
  "cellMG-GMM"             = '"cellMG-GMM"',
  "cellGMM"                = "cellGMM"
)

g_ell = ggplot() +
    geom_polygon(data = ellipse_data,
                 aes(x = get(var1),
                     y = get(var2),
                     group =  interaction(method, group),
                     col = as.factor(group)
                     ),
                 fill = "transparent",
                 alpha = 0.05) +
  geom_point(data = data_plot,
             aes(x = get(var1),
                 y = get(var2),
                 col = as.factor(class)),
             alpha = 0.05,
             shape = 16) +
  facet_grid(rows = vars(method),
             cols = vars(group),
             labeller = as_labeller(quality_labels, label_parsed)) +
  theme_bw(base_size = 16) +
  labs(x = "alcohol", y = "residual sugar") +
  scale_color_manual(name = "", values = cols, labels = c(
    "1" = "Low Quality (Model)",
    "2" = "Middle Quality (Model)",
    "3" = "High Quality (Model)")) +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0))
g_ell
ggsave(g_ell, file = "wine_white_ellipses.pdf", width = 10, height = 8)


