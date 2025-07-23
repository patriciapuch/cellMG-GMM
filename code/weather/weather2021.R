library(ssMRCD)
library(dplyr)
library(ggplot2)
library(stringr)
library(rnaturalearth)
library(tidyr)
library(rnaturalearthdata)
library(ggnewscale)
library(ellipse)
library(geomtextpath)
library(ggrepel)


# prepare data from ssMRCD package
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


# construct groups based on spatial proximity
cut_lon = c(min(weatherAUT2021$lon)-0.2, 12, 16, max(weatherAUT2021$lon)+ 0.2)
cut_lat = c(min(weatherAUT2021$lat)-0.2, 48, max(weatherAUT2021$lat)+0.2)
groups = ssMRCD::groups_gridbased(weatherAUT2021$lon, 
                                  weatherAUT2021$lat, 
                                  cut_lon, cut_lat)
table(groups)
N = length(unique(groups))


# calculate model
start = Sys.time()
out = ssMRCD::cellMGGMM(X = data,
                        groups = groups,
                        nsteps = 100,
                        alpha = 0.5,
                        hperc = 0.75,
                        maxcond = 100, 
                        silent = TRUE)
end = Sys.time()
end-start


# load Austria as sf object
austria <- ne_countries(scale = "medium", country = "Austria", returnclass = "sf")

g_boundary = ggplot() + 
  geom_sf(data = austria, fill = "transparent", color = "black") +
  theme_classic()

g_mapgrid = g_boundary +
  geom_point(aes(x = weatherAUT2021$lon, 
                 y = weatherAUT2021$lat,
                 shape = as.factor(apply(out$probs,1, which.max))), 
             col = "black") +
  geom_hline(aes(yintercept = cut_lat), 
             linetype = 2) + 
  geom_vline(aes(xintercept = cut_lon),
             linetype = 2) + 
  labs(title = "",
       color = "m",
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


# results
out$rho
plot(out$objvals)
# mixing probabilities
cat("Pi (in %):\n")
round(out$pi_groups*100, 2)


# percentage of outliers
cat("% Outliers per group and variable:\n")
round(sapply(1:N, function(x) colMeans(1-out$W[groups == x, ]))*100, 2)


# calculate residuals
res = residuals_mggmm(X = data, 
                      groups = groups,
                      Sigma = out$Sigma,
                      mu = out$mu, 
                      probs = out$probs,
                      W = out$W, 
                      set_to_zero = TRUE)


# make residual map with probabilities
subset_resids = tidyr::pivot_longer(cbind(data.frame(res), 
                                          stations_short, 
                                          groups = as.character(groups),
                                          probs = out$probs)[rowSums(1-out$W) > 0 , ],
                                    cols = c(p:rel, probs.1:probs.5)) %>%
  mutate(name = gsub("probs.", "Group ", name),
         groups = paste0("Group ", groups),
         groups_probs = ifelse(grepl("Group", name), groups, "vv"),
         value_probs = ifelse(grepl("Group", name), value, NA),
         value_resids = ifelse(!grepl("Group", name), value, NA),
         facet =  ifelse(grepl("Group", name), "Probabilities", "Residuals"))


# plot outlying cells
g_cell = ggplot(subset_resids) + 
  geom_tile(aes(y = stations,
                x = name, 
                fill = value_resids), 
            col = "white") + 
  scale_fill_gradient2("Resid.", 
                       low = scales::muted("blue"), 
                       high = scales::muted("red"),
                       na.value = "transparent") +
  new_scale_fill()+
  geom_tile(aes(y = stations,
                x = name, 
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


# bivariate ellipse plot
var1 = "t"
var2 = "vv"
which.var = which(colnames(data) %in% c(var1, var2))

data_plot = cbind(data, stations_short)[rowSums(1-out$W[, which.var]) != 0, ]

data_shape = rep(NA, n)
data_shape[rowSums(out$W[, which.var]) == 2] = "No outliers"
data_shape[out$W[, which.var[1]] == 0] = paste0("Outliers in ", var1)
data_shape[out$W[, which.var[2]] == 0] = paste0("Outliers in ", var2)
data_shape[rowSums(1-out$W[, which.var]) == 2] = "Outliers in both"
data_shape = data_shape[rowSums(1-out$W[, which.var]) != 0]


# plot ellipses
g_ells = ggplot()
lty = c("solid", "aa", "dashed",  "dotted", "dotdash")
h = c(0.5, 0.5, 0.5, 0.4, 0.85)
for(i in 1:N){
  cov = out$Sigma[[i]][c(var1, var2), c(var1, var2)]
  mean = out$mu[[i]][c(var1, var2)]
  ellipse_data <- as.data.frame(ellipse(cov, centre = mean, level = 0.95))
  
  
  # cov = out_smoothed$Sigma[[i]][c(var1, var2), c(var1, var2)]
  # mean = out_smoothed$mu[[i]][c(var1, var2)]
  # ellipse_data <- as.data.frame(ellipse(cov, centre = mean, level = 0.99))
  
  g_ells = g_ells + 
    geom_textpath(data = ellipse_data, 
                  aes(x = get(var1), 
                      y = get(var2)),
                  label = paste("Group", i),
                  size = 4,
                  alpha = 0.5,
                  hjust = h[i]
    ) +
    geom_polygon(data = ellipse_data,
                 aes(x = get(var1),
                     y = get(var2)),
                 fill = scales::alpha("black", alpha = 0.05))
  
}
g_ells


# add observations
set.seed(1)
g_bivariate = g_ells +
  geom_point(data = data_plot,
             aes(x = get(var1),
                 y = get(var2),
                 #color = as.factor(groups[rowSums(1-out$W[, which.var]) != 0]),
                 shape = as.factor(groups[rowSums(1-out$W[, which.var]) != 0])
             ),
             #alpha = 0.5,
             size = 2,
             stroke = 1) +
  #scale_color_manual(name = "Group", values = colors) +
  scale_shape_manual(name = "Group", values = c(16, 17, 15, 3, 4)) +
  new_scale_color() +
  geom_label_repel(data = data_plot,
                   aes(x = get(var1),
                       y = get(var2),
                       label = stations,
                       col = data_shape),
                   size = 4,
                   fill = "white" #scales::alpha(colour = "white", alpha = 0.5)
                   #max.overlaps = 5,
  )+
  geom_point(data = data_plot,
             aes(x = get(var1),
                 y = get(var2),
                 col = data_shape),
             alpha = 0
  )+
  scale_color_manual("Outlying in", 
                     values = c("mediumpurple", 
                                scales::muted("red"), 
                                scales::muted("cornflowerblue")),
                     labels = c("t and vv", "vv", "t")) + 
  guides(colour = guide_legend(override.aes = aes(label = "")))+
  theme_classic(base_size = 14)+
  xlim(c(-7, 14)) +
  ylim(c(0, 8.5)) +
  #scale_y_log10()+
  labs(x = var1, y = var2, title = "") +
  theme(aspect.ratio = 1/2.5)

g_bivariate