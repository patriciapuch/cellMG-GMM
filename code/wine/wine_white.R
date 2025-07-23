# WHITE WINE DATA (UCI Machine Learning Repository)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ssMRCD)
library(cellWise)

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

# apply cellMG-GMM method to data
out = ssMRCD::cellMGGMM(X = wine_white,
                       groups = groups_white,
                       nsteps = 100,
                       alpha = 0.75,
                       maxcond = 100)

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

# parallel coordinate plot - preparation
dat_all = rbind(cbind(data.frame(wine_white),
                      kind = "data",
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
colnames(out$W) = paste0("W", colnames(out$W))
stand_d_sub = dat_all %>% filter(kind == "data")
stand_d_sub[, 1:11][out$W == 1] = NA


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
                alpha = interaction(assign, groups)),
            linewidth = 0.5) +
  geom_line(data = datplot %>% filter(kind == "mu"),
            aes(x = name,
                y = value,
                #col = as.factor(groups),
                group = id),
            col = "black",
            linewidth = 0.75) +
  geom_errorbar(data = datplot %>%
                filter(kind != "data") %>%
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
#ggsave(g_par, file = "wine_white_pcp.pdf", width = 10, height = 6)
