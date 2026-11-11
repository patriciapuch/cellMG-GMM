# Weather Example (Appendix)
# set your working directory in the same folder as weather.R

library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)
library(grid) # for arrows
library(plotly)
library(cellWise)
library(rrcov)
library(ssMRCD)
library(sf)
library(ggnewscale)
library(stringr)
library(ellipse)
library(patchwork)
library(ggnewscale)
library(ggrepel)
library(geomtextpath)

source("../simulations/sim_setup.R")
set.seed(123)

# prepare data
data("weatherAUT2021")

data = weatherAUT2021 %>% select(p:rel)
stations = weatherAUT2021$name
n = dim(data)[1]

stations_short = data.frame(stations) %>%
  mutate(stations = case_when(stations == "MARIAZELL/ST.SEBASTIAN-FLUGFELD" ~ "MARIAZELL",
                              stations == "RAX/SEILBAHNBERGSTATION" ~ "RAX/BERGST",
                              stations == "HOHE WAND/HOCHKOGELHAUS" ~ "HOHE WAND",
                              stations == "SONNBLICK (TAWES)" ~ "SONNBLICK",
                              stations == "ACHENKIRCH CAMPINGPLATZ" ~ "ACHENKIRCH",
                              stations == "LINZ-STADT" ~ "LINZ",
                              stations == "WIEN-INNERE STADT" ~ "WIEN-I.",
                              stations == "WIEN-JUBILAEUMSWARTE" ~ "WIEN-J.",
                              stations == "JAUERLING/ORF" ~ "JAUERLING",
                              TRUE ~ stations)) %>%
  mutate(stations = str_to_title(tolower(stations))) %>%
  mutate(stations = factor(stations, levels = sort(stations, decreasing = TRUE), ordered = TRUE ))

## group by longitude-latitude (N = 5)
cut_lon = c(min(weatherAUT2021$lon)-0.2, 12, 16, max(weatherAUT2021$lon)+ 0.2)
cut_lat = c(min(weatherAUT2021$lat)-0.2, 48, max(weatherAUT2021$lat)+0.2)
groups = ssMRCD::groups_gridbased(weatherAUT2021$lon,
                                  weatherAUT2021$lat,
                                  cut_lon, cut_lat)
N = length(unique(groups))

# hyper-parameter tuning for alpha
alpha_grid = seq(0.5, 1, 0.05)
l = vector(mode = "list", length = length(alpha_grid))
for(i in 1:length(alpha_grid)){
  l[[i]] = ssMRCD::cellMGGMM(X = data,
                          groups = groups,
                          nsteps = 100,
                          alpha = alpha_grid[i],
                          maxcond = 100)
}
out_tune = select_alpha(list_results = l, alpha_grid = alpha_grid)
obj = out_tune[-c(1,2)]

# plot hyperparameter tuning
g_alpha = ggplot() +
  geom_line(aes(x = alpha_grid, y = obj)) +
  geom_point(aes(x = out_tune[1], y = obj[out_tune[2]]), col = "red") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust= 1, vjust = 0.5),
        panel.grid = element_blank(),
        aspect.ratio = 1/2) +
  labs(x = expression(alpha),
       y = "Objective function")
ggsave(g_alpha, file = "weather_alphatuning.pdf", width = 5)


# calculate model
out = ssMRCD::cellMGGMM(X = data,
               groups = groups,
               nsteps = 100,
               alpha = out_tune[1],  # 0.7
               hperc = 0.75,
               maxcond = 100)

# pi matrix
cat("Pi (in %):\n")
round(out$pi_groups*100, 2)

# percentage of outliers
cat("% Outliers per group and variable:\n")
round(sapply(1:N, function(x) colMeans(1-out$W[groups == x, ]))*100, 2)

# plot map of weather stations
load("Hoehenprofil_AUT_Dataframe.RData")
boundary_BL <- st_read("Bundeslaender_50.shp")

colors = c("#283618", "#606C38", "#FEFAE0", "#DDa15E", "#BC6C25")
g_mapgrid = ggplot() +
  geom_sf(data = boundary_BL,
          alpha = 0,
          fill = "transparent",
          col = "transparent") +
  geom_raster(data = grey_proj_dat,
              aes(x = x,
                  y = y,
                  fill = (dhm_at_lamb_10m_2018)),
              alpha = 0.7) +
  geom_point(aes(x = weatherAUT2021$lon,
                 y = weatherAUT2021$lat,
                 shape = as.factor(apply(out$probs,1, which.max))),
             col = "black") +
  geom_hline(aes(yintercept = cut_lat),
             linetype = 2) +
  geom_vline(aes(xintercept = cut_lon),
             linetype = 2) +
  labs(title = "",
       shape = "Groups",
       x = "",
       y = "") +
  theme_classic(base_size = 13) +
  scale_shape_manual(values = c(16, 17, 15, 3, 4)) +
  scale_fill_gradientn(colors = colors,
                       breaks = c(1000, 2000, 3000),
                       labels = c(1,2,3),
                       values = c(0, 0.08, 0.12, 0.24, 0.366, 0.55, 1 ),
                       name = "Altitude \n[1000m]", ) +
  theme(legend.position = "bottom")
g_mapgrid
ggsave(g_mapgrid, file = "weather_map.pdf", height = 5)



# calculate residuals
res = residuals_mggmm(X = data,
                     groups = groups,
                     Sigma = out$Sigma,
                     mu = out$mu,
                     probs = out$probs,
                     W = out$W,
                     set_to_zero = TRUE)


# make residual - posterior probability plot
subset_resids = melt(cbind(data.frame(res),
                           stations_short,
                           groups = as.character(groups),
                           probs = out$probs)[rowSums(1-out$W) > 0 , ],
                     id.vars = c("stations", "groups"))  %>%
  mutate(variable = gsub("probs.", "Group ", variable),
         groups = paste0("Group ", groups),
         groups_probs = ifelse(grepl("Group", variable), groups, "vv"),
         value_probs = ifelse(grepl("Group", variable), value, NA),
         value_resids = ifelse(!grepl("Group", variable), value, NA),
         facet =  ifelse(grepl("Group", variable), "Probabilities", "Residuals"))

# plot outlying cells
g_cell = ggplot(subset_resids) +
  geom_tile(aes(y = stations,
                x = variable,
                fill = value_resids),
                col = "white") +
  scale_fill_gradient2("Resid.",
                       low = scales::muted("blue"),
                       high = scales::muted("red"),
                       na.value = "transparent") +
  new_scale_fill()+
  geom_tile(aes(y = stations,
                x = variable,
                fill = value_probs),
            col = "white") +
  scale_fill_gradient2("Prob.",
                       low = "white",
                       high = scales::muted("blue4"),
                       midpoint = 0.5,
                       mid = scales::muted("cornflowerblue"),
                       na.value = "transparent") +
  geom_point(aes(y = stations,
                 x = groups_probs,
                 col = facet)) +
  scale_color_manual("",values = c("grey", "transparent")) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust= 1, vjust = 0.5),
        panel.grid = element_blank(),
        aspect.ratio = 3.7) +
  labs(x = "",
       y = "") +
  facet_grid(cols = vars(facet),
             scales = "free") +
  guides(col = "none")
g_cell
ggsave(g_cell, file = "weather_residuals.pdf", height = 6)


# make ellipse plot
var1 = "t"
var2 = "vv"
which.var = c(which(colnames(data) == var1), which(colnames(data) == var2))
data_plot = cbind(data, stations_short)[rowSums(1-out$W[, which.var]) != 0, ]



data_shape = rep(NA, n)
data_shape[rowSums(out$W[, which.var]) == 2] = "No outliers"
data_shape[out$W[, which.var[1]] == 0] = "Outliers in t"
data_shape[out$W[, which.var[2]] == 0] = "Outliers in vv"
data_shape[rowSums(1-out$W[, which.var]) == 2] = "Outliers in both"
data_shape = data_shape[rowSums(1-out$W[, which.var]) != 0]


g_ells = ggplot()
h = c(0.5, 0.5, 0.5, 0.4, 0.85)
for(i in 1:N){
  cov = out$Sigma[[i]][c("t", "vv"), c("t", "vv")]
  mean = out$mu[[i]][c("t", "vv")]
  ellipse_data <- as.data.frame(ellipse(cov, centre = mean, level = 0.95))

  g_ells = g_ells +
    geom_textpath(data = ellipse_data,
            aes(x = t,
                y = vv),
            label = paste("Group", i),
            size = 4,
            alpha = 0.5,
            hjust = h[i]
            ) +
    geom_polygon(data = ellipse_data,
            aes(x = t,
                y = vv),
            fill = scales::alpha("black", alpha = 0.05))
}

set.seed(1)
g_bivariate = g_ells +
  geom_point(data = data_plot,
             aes(x = t,
                y = vv,
                shape = as.factor(groups[rowSums(1-out$W[, which.var]) != 0])
                ),
             size = 2,
             stroke = 1) +
  scale_shape_manual(name = "Group", values = c(16, 17, 15, 3, 4)) +
  new_scale_color() +
  geom_label_repel(data = data_plot,
            aes(x = t,
                y = vv,
                label = stations,
                col = data_shape),
            size = 4,
            fill = "white"
            )+
  geom_point(data = data_plot,
            aes(x = t,
                y = vv,
                col = data_shape),
            alpha = 0
            )+
  scale_color_manual("Outlying in",
                     values = c("mediumpurple",
                                scales::muted("red"),
                                scales::muted("cornflowerblue")),
                     labels = c("t and vv", "t", "vv"))+
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  theme_classic(base_size = 14)+
  xlim(c(-7, 14)) +
  ylim(c(0, 8.5)) +
  labs(x = var1, y = var2, title = "") +
  theme(aspect.ratio = 1/2.5)
g_bivariate
ggsave(g_bivariate, file = "weather_tempvswind.pdf", width = 10)
