library(robustmatrix)
library(dplyr)
library(ggplot2)
library(cellWise)
library(reshape2)
library(tidyr)
library(ggh4x)
library(patchwork)
library(ssMRCD)

# load data from robustmatrix package
data(darwin)

# construct features
darwin_avg = apply(darwin, c(3,1), median)
darwin_avg = darwin_avg[,! colnames (darwin_avg) %in% c("air_time", "total_time", "mean_gmrt")]
colnames(darwin_avg) = paste0(colnames(darwin_avg) , "_med")
colnames(darwin_avg) = colnames(darwin_avg) %>% gsub("_in_", "_", .) %>% gsub("_on_", "_", .)

darwin_mad = apply(darwin, c(3,1), mad)
darwin_mad = darwin_mad[,! colnames (darwin_mad) %in% c("air_time", "total_time", "mean_gmrt")]
colnames(darwin_mad) = paste0(colnames(darwin_mad) , "_mad")
colnames(darwin_mad) = colnames(darwin_mad) %>% gsub("_in_", "_", .) %>% gsub("_on_", "_", .)

darwin_all = cbind(darwin_avg, darwin_mad)

# scale data for residual calculation
scale_loc = cellWise::estLocScale(darwin_all, type = "mcd")
darwin_sca = scale(darwin_all, center = scale_loc$loc, scale = scale_loc$scale)

# construct groups
groups = as.numeric(as.factor(rownames(darwin_avg)))
N = max(groups)

# run method for varying alpha levels
lambda_seq = seq(0.5, 1, 0.01)
out_all = list()
probs = matrix(NA, ncol = length(groups), nrow = length(lambda_seq))
rownames(probs) = lambda_seq

resids = array(NA,
               dim = c(length(groups), length(lambda_seq), dim(darwin_sca)[2]),
               dimnames = list(1:174, lambda_seq, colnames(darwin_sca)))

outs_frac1 = matrix(NA, ncol = dim(darwin_sca)[2], nrow = length(lambda_seq))
rownames(outs_frac1) = lambda_seq

outs_frac2 = matrix(NA, ncol = dim(darwin_sca)[2], nrow = length(lambda_seq))
rownames(outs_frac2) = lambda_seq

for(i in 1:length(lambda_seq)){
  out = ssMRCD::cellMGGMM(X = darwin_sca,
                         groups = groups,
                         nsteps = 100,
                         alpha = lambda_seq[i],
                         maxcond = 100)
  out_all = c(out_all, list(out))
  probs[i, groups == 1] = out$probs[groups == 1, 1]
  probs[i, groups == 2] = out$probs[groups == 2, 2]

  outs_frac1[i, ] = colSums(1-out$W[groups == 1, ])/sum(groups == 1)
  outs_frac2[i, ] = colSums(1-out$W[groups == 2, ])/sum(groups == 2)

  res = residuals_mggmm(X = darwin_sca,
                       groups = groups,
                       Sigma = out$Sigma,
                       mu = out$mu,
                       probs = out$probs,
                       W = out$W,
                       set_to_zero = TRUE)
  resids[, i, ] = res
}

# probability plots (left panel)
id_change = which(colSums((1-probs) > 0.01) != 0)

tmp = apply(probs, 2, function(x)  max(which( x < 0.8)))
sorted_ids = sort.int(tmp, index.return = TRUE, decreasing = TRUE)$ix

col = 1:174
col[!(1:174 %in% id_change)] = NA

plot_probs = cbind(data.frame(t(probs)),
                   groups = ifelse(groups == 1, "Healthy", "Alzheimer"),
                   id = 1:174,
                   col = col) %>%
  data.frame %>%
  tidyr::pivot_longer(X0.5:X1) %>%
  mutate(alpha = as.numeric(gsub("X", "", name)))

g_probs1 = ggplot() +
  geom_tile(data = plot_probs[plot_probs$id %in% id_change & plot_probs$groups == "Alzheimer", ],
            aes(x = alpha,
                fill = value,
                col = value,
                y = factor(id, levels = sorted_ids))) +
  scale_x_reverse(expand = c(0, 0), breaks = c(1, 0.75, 0.5)) +
  scale_y_discrete(expand = c(0, 0.85))+
  facet_wrap(vars(groups))+
  theme_classic(base_size = 13) +
  scale_fill_gradient(name = expression(t[gig]),
                      low = "grey80", high = "grey30",
                      breaks = c(0, 0.5, 1),
                      labels = c(0, 0.5, 1)) +
  scale_color_gradient(name = expression(t[gig]),
                       low = "grey80", high = "grey30",
                       breaks = c(0, 0.5, 1),
                       labels = c(0, 0.5, 1)) +
  labs(y = "", x = expression(alpha)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        aspect.ratio = 2*17/10,
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))

g_probs2 = ggplot() +
  geom_tile(data = plot_probs[plot_probs$id %in% id_change & plot_probs$groups == "Healthy", ],
            aes(x = alpha,
                fill = value,
                col = value,
                y = factor(id, levels = sorted_ids))) +
  scale_x_reverse(expand = c(0, 0), breaks = c(1, 0.75, 0.5)) +
  scale_y_discrete(expand = c(0, 0.85)) +
  facet_wrap(vars(groups)) +
  theme_classic(base_size = 13) +
  scale_fill_gradient(name = expression(t[gig]),
                      low = "grey80", high = "grey30",
                      breaks = c(0, 0.5, 1),
                      labels = c(0, 0.5, 1)) +
  scale_color_gradient(name = expression(t[gig]),
                       low = "grey80", high = "grey30",
                       breaks = c(0, 0.5, 1),
                       labels = c(0, 0.5, 1)) +
  labs(y = "", x = expression(alpha)) +
  theme(legend.position = "bottom",
        aspect.ratio = 2*10/10,
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  guides(color = guide_colorbar(barwidth  = unit(2, "cm"))) +
  guides(fill = guide_colorbar(barwidth  = unit(2, "cm")))
g_probs1 / g_probs2


# residual variation plot
ind_noouts = apply(abs(resids), c(1,3), sum) == 0
ind_noouts_sum = apply(resids != 0, c(1,3), sum) %>%
  cut( breaks = c(0, 5, 51), include.lowest = TRUE) %>%
  matrix(ncol = 30)

resids_sd = apply(resids, c(1,3), sd) %>%
  data.frame %>%
  mutate(groups = ifelse(groups == 1, "Healthy", "Alzheimer"),
         id = factor(1:174, level = sorted_ids))
resids_sd[, 1:30][ind_noouts] = NA
resids_sd = resids_sd %>%
  melt(id.vars = c("groups", "id")) %>%
  mutate(change_id = factor(ifelse(id %in% id_change, "Switcher", "Stable"),
                            level = c("Switcher", "Stable")),
         variable = factor(variable,
                           level = paste0(rep(gsub("_med", "",colnames(darwin_avg)), each = 2),
                                          rep(c("_med", "_mad"), times = 15))))


resids_sd = cbind(resids_sd, nouts = melt(ind_noouts_sum)$value)
ind_allouts = which(apply(resids != 0, c(1), sum) != 0)  # all cells not outlying (84 in plot, 90 are not shown)
rowSums(1-ind_noouts)


g_residsd = ggplot(resids_sd %>% filter(id %in% c(ind_allouts, id_change))) +
  geom_tile(aes(x = variable,
                y = id,
                fill = value),
            col = "grey90") +
  geom_point(aes(x = variable,
                y = id,
                shape = as.factor(nouts),
                col = as.factor(nouts)),
             stroke = 0.9,
             size = 1) +
  scale_fill_gradient(name = "Residual\nvariation",
                      na.value = "white",
                      low ="lightblue1",
                      high = "lightblue4",
                      limits = c(0, 2.2),
                      breaks = c(0, 0.5, 1, 1.5, 2),
                      labels = c(0, 0.5, 1, 1.5, 2)) +
  scale_shape_manual(name = "Number of\noutlying\ncells over \U03B1",
                     values = c(19, NA),
                     label = c("[6, 51]", "[0,5]")) +
  scale_color_manual(name = "Number outlying \ncells over \U03B1",
                     values = c ("black","black",3), guide = "none") +
  theme_classic(base_size = 13)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1),
        axis.text.y = element_text(size = 9),
        legend.position = "right",plot.margin = margin(0, 0, 0, 0)) +
  facet_nested(groups + change_id ~.,
               scales = "free",
               space = "free"
               )+
  labs(y = "ID",
       x = "")
g_residsd
