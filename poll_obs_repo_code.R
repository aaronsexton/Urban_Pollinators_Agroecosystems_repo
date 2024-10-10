
# Repository code for analyses conducted relating to 
# "Urban pollinator communities are structured by local-scale garden features, not landscape context"

# Code originally archived on Oct.10th, 2024
#  by Aaron N. Sexton



# Pollinator Community Composition ----

library(vegan)
library(tidyverse)


# First we will run the NMDS and vector analysis on pollinator community composition

# Load in pollinator observations
poll_obs <- read.csv("../Data/Poll_Obs_repos/pollinator_observations.csv")
str(poll_obs)

# Join the observations with the features data
features <- read.csv("../Data/Poll_Obs_repos/garden_features.csv")
str(features)

# Join them
poll_obs <- left_join(features, poll_obs)
str(poll_obs)

# Build a distance matrix of the pollinator values
poll_dist <- vegdist(decostand(poll_obs[17:36], "norm"), "euclidean")

# Run the NMDS
poll_nms <- metaMDS(poll_dist, distance = "euclidean", k = 3, maxit = 999,
                    trymax = 250, wascores = T)
poll_nms
plot(poll_nms, "sites")   # Produces distance 
orditorp(poll_nms, "sites")   # Gives points labels
stressplot(poll_nms)

# Create a predictor df to run the vector analysis
preds      <- poll_obs[6:16]
preds$city <- poll_obs$city

# Run the vector analysis
en <- envfit(poll_nms, preds, permutations = 999, na.rm = TRUE)
en
plot(poll_nms)
plot(en)

# Save the coordinates of sites and species and vectors for nicer plotting
data.scores <- as.data.frame(poll_nms$points)
en_coord <- as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
poll_nms$species
poll_nms$nobj
spp.rel <- decostand(poll_obs[17:36], "norm")
spp_scrs <- 
  sppscores(poll_nms) <- spp.rel


plot(poll_nms)
plot(poll_nms, type = "text")

nms_species <- as.data.frame(poll_nms$species)
nms_species$Taxa <- rownames(nms_species)

data.scores$city <- poll_obs$city


# Make names pretty
en_coord$vars  <- rownames(en_coord)
en_coord[1,3]  <- "Veg Height"
en_coord[2,3]  <- "Floral Richness"
en_coord[3,3]  <- "Sand Cover"
en_coord[4,3]  <- "Bare Soil" 
en_coord[5,3]  <- "Floral Abundance"
en_coord[6,3]  <- "C/W Ratio"
en_coord[7,3]  <- "Woody Abundance"
en_coord[8,3]  <- "Woody Richness"
en_coord[9,3]  <- "Canopy Cover"
en_coord[10,3] <- "Garden Size"
en_coord[11,3] <- "Urbanization"

en_coord$Significance <- "ns"
en_coord[1,4]  <- "< 0.01"
en_coord[2,4]  <- "< 0.01"
en_coord[4,4]  <- "< 0.01"
en_coord[8,4]  <- "< 0.01"


# Plot out the nmds for both the continuous vectors and the categorical vector (city) 


city_nmds <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) +
  geom_segment(data = en_coord, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size = 1, alpha = 0.001) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("cornflowerblue", "tomato3")) +
  geom_label(data = nms_species, aes(x = MDS1, y = MDS2, label = Taxa),
             fontface = "bold", alpha = 0.8, size = 4, color = "black") + 
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        legend.key = element_blank(), 
        legend.title = element_text(size = 24, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 22, colour = "grey30"))
city_nmds


vector_nmds <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(data = data.scores, aes(size = MDS3), alpha = 0.6, color = "mediumpurple1") + 
  geom_segment(data = en_coord, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2,
                                     color = Significance), size = 1) +
  geom_text(data = en_coord, aes(x = (NMDS1 + 0.02), y = (NMDS2 - 0.01),
                                  color = Significance), 
            fontface = "bold", label = en_coord$vars) + 
  scale_color_manual(values = c("black", "gray65")) +
  geom_label(data = nms_species, aes(x = MDS1, y = MDS2, label = Taxa),
             fontface = "bold", alpha = 0.8, size = 4, color = "black") + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        legend.key = element_blank(), 
        legend.title = element_text(size = 24, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 22, colour = "grey30"))
vector_nmds



ggpubr::ggarrange(vector_nmds, city_nmds)




rm(list = ls())







# Calculate Beta diversity ----

library(vegan)
library(tidyverse)



# Read in the dataframe with all pollinator observation
poll_obs_fg <- read.csv("../Data/Poll_Obs_repos/pollinator_observations.csv")
str(poll_obs_fg)

# Create site and city objects for beta diversity calculations
sites <- poll_obs_fg$garden_code
city <- poll_obs_fg$city

# Run distance matrix for pollinator observations
poll_dist <- vegdist(decostand(poll_obs_fg[6:25], "norm"), "euclidean")

# Calculate beta diversity at the site level
bd <- betadisper(poll_dist, sites)
plot(bd)
bd$group.distances
hist(bd$group.distances)

# Save values in a dataframe
site_betas <- as.data.frame(bd$group.distances)
site_betas$garden_code <- rownames(site_betas)
site_betas <- rename(site_betas, pollinator_beta = `bd$group.distances`)

# Run it also at the city level
bd_city <- betadisper(poll_dist, city)
plot(bd_city)
bd_city$group.distances
anova(bd_city)
# Berlin    Munich 
# 0.6861744 0.6778099 
# No difference in dispersion between cities



# Run regressions ----

# Read in a pre-existing database with plant beta and alpha diversity values
floral_betas <- read.csv("../Data/Poll_Obs_repos/floral_divers.csv")
str(floral_betas)

# Combine beta div files
betas2 <- left_join(floral_betas, site_betas)



# Pollinator beta ~ Urbanization
a <- ggplot(betas2, aes(x = impervious_1000, y = pollinator_beta)) + 
  geom_point(aes(color = city), size = 5, alpha = 0.8) + 
  theme_bw() +
  labs(color = "City", x = "Impervious Surface (1000m radius)", y = "Pollinator Beta Diversity") +
  scale_color_manual(values = c("cornflowerblue", "tomato3")) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 32, hjust = 0.5)) +
  annotate("text", x = 12, y = 0.79, label = "a)", size = 12)
a 

poll_urb_lm <- lm(pollinator_beta ~ scale(impervious_1000), 
                  data = betas2)
summary(poll_urb_lm)


# Pollinator Beta  ~ Floral Beta
b <- ggplot(betas2, aes(x = floral_beta, y = pollinator_beta)) +
  geom_point(aes(color = city), size = 5, alpha = 0.8) +
  theme_bw() +
  labs(color = "City", x = "Floral Beta Diversity", y = "Pollinator Beta Diversity") +
  scale_color_manual(values = c("cornflowerblue", "tomato3")) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 32, hjust = 0.5)) +
  annotate("text", x = 0.765, y = 0.79, label = "c)", size = 12)
b

summary(lm(betas2$pollinator_beta ~ betas2$floral_beta))



# Floral Beta ~ Urbanization
fb_u_lm <- ggplot(betas2, aes(x = impervious_1000, y = floral_beta)) + 
  geom_smooth(method = "lm", alpha = 0.1, color = "purple3", fill = "purple4") +
  geom_point(aes(color = city), size = 5, alpha = 0.8) + 
  theme_bw() +
  labs(color = "City", x = "Impervious Surface (1000m radius)", y = "Floral Beta Diversity") +
  scale_color_manual(values = c("cornflowerblue", "tomato3")) +
  scale_fill_manual(values = c("cornflowerblue", "tomato3")) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 32, hjust = 0.5)) +
  annotate("text", x = 12, y = 0.95, label = "b)", size = 12)
fb_u_lm

flB_urb_lm2 <- lm(floral_beta ~ scale(impervious_1000), 
                  data = betas2)
summary(flB_urb_lm2)


# Pollinator beta ~ Floral alpha
pb_fa_lm <- ggplot(betas2, aes(x = plant_rich, y = pollinator_beta)) + 
  geom_smooth(method = "lm", alpha = 0.1, color = "purple3", fill = "purple4") +
  geom_point(aes(color = city), size = 5, alpha = 0.8) + 
  theme_bw() +
  labs(color = "City", x = "Floral Alpha Richness", y = "Pollinator Beta Diversity") +
  scale_color_manual(values = c("cornflowerblue", "tomato3")) +
  scale_fill_manual(values = c("cornflowerblue", "tomato3")) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 32, hjust = 0.5)) +
  annotate("text", x = 3, y = 0.79, label = "d)", size = 12)
pb_fa_lm

summary(lm(betas2$pollinator_beta ~ betas2$plant_rich))




# Put them all together
ggpubr::ggarrange(a, fb_u_lm, b, pb_fa_lm,
                  ncol = 2, nrow = 2, 
                  common.legend = T,
                  legend = "right")










# City-level comparisons ----


# Need to calculate average pollinator alpha div at the site level

# Group pollinator observation data by site
poll_obs_grpd_plot <- dplyr::group_by(poll_obs_fg, garden_code)
str(poll_obs_grpd_plot)

# Calculate average value
sum1 <- dplyr::summarise(poll_obs_grpd_plot, 
                         poll_rich = mean(richness_poll))
hist(sum1$poll_rich)

# Add it to the Beta diversity file
betas2 <- left_join(betas2, sum1)


# Check distributions of alpha and beta diversity
hist(betas2$poll_rich)
hist(betas2$plant_rich)
hist(betas2$floral_beta)
hist(betas2$pollinator_beta)

# Need to run a Shaprio-Wilk test to check for normality
# First split by city then test each
mun_betas <- filter(betas2, city == "Munich")
ber_betas <- filter(betas2, city == "Berlin")

# Munich tests
shapiro.test(mun_betas$plant_rich) # Normal
shapiro.test(mun_betas$poll_rich) # Normal
shapiro.test(mun_betas$floral_beta)  # Not normal
shapiro.test(mun_betas$pollinator_beta) # Normal
# Berlin tests
shapiro.test(ber_betas$plant_rich) # Normal
shapiro.test(ber_betas$poll_rich) # Normal
shapiro.test(ber_betas$floral_beta)  # Normal
shapiro.test(ber_betas$pollinator_beta) # Normal

# Okay so the only non-normal is Munich Floral Beta div

# Run tests

# Kruskal-Wallis test for floral beta
kruskal.test(floral_beta ~ city, data = betas2)
# Non-significant
# Chi-squared = 0.219, p = 0.639


# T-test for the rest

# Plant Alpha
plant_alpha_t <- t.test(plant_rich ~ city, data = betas2)
plant_alpha_t
# t = 0.27768, p-value = 0.7831

# Pollinator Alpha
poll_alpha_t <- t.test(poll_rich ~ city, data = betas2)
poll_alpha_t
# t = 1.4405, p-value = 0.162

# Pollinator Beta
poll_beta_t <- t.test(pollinator_beta ~ city, data = betas2)
poll_beta_t
# t = 0.5523, p-value = 0.5858






## Violin plots ----



a <- ggplot(betas2, aes(x = city, y = pollinator_beta)) +
  geom_violin(aes(fill = city), trim = F, alpha = 0.8) +
  geom_jitter(position = position_jitter(0.1), alpha = 0.2, size = 4) +
  theme_bw() +
  labs(y = "Pollinator Beta Diversity") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 24)) +
  scale_fill_manual(values = c("cornflowerblue", "tomato3")) +
  annotate("text", x = 0.5, y = 0.95, label = "a)", size = 12)
a

b <- ggplot(betas2, aes(x = city, y = floral_beta)) +
  geom_violin(aes(fill = city), trim = F, alpha = 0.8) +
  geom_jitter(position = position_jitter(0.1), alpha = 0.2, size = 4) +
  theme_bw() +
  labs(y = "Floral Beta Diversity") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 24)) +
  scale_fill_manual(values = c("cornflowerblue", "tomato3")) +
  annotate("text", x = 0.5, y = 1, label = "b)", size = 12)
b

d <- ggplot(betas2, aes(x = city, y = plant_rich)) +
  geom_violin(aes(fill = city), trim = F, alpha = 0.8) +
  geom_jitter(position = position_jitter(0.1), alpha = 0.2, size = 4) +
  theme_bw() +
  labs(y = "Floral Alpha Diversity") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 24)) +
  scale_fill_manual(values = c("cornflowerblue", "tomato3")) +
  annotate("text", x = 0.5, y = 21, label = "d)", size = 12)


c <- ggplot(betas2, aes(x = city, y = poll_rich)) +
  geom_violin(aes(fill = city), trim = F, alpha = 0.8) +
  geom_jitter(position = position_jitter(0.1), alpha = 0.2, size = 4) +
  theme_bw() +
  labs(y = "Pollinator Alpha Diversity") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 24)) +
  scale_fill_manual(values = c("cornflowerblue", "tomato3")) +
  annotate("text", x = 0.5, y = 7.2, label = "c)", size = 12)
c




ggpubr::ggarrange(a, b, c, d,
                  ncol = 2, nrow = 2)





rm(list = ls())








# Mapping ----

library(sf)
library(terra)
library(tidyverse)
library(ggspatial)


#### Berlin ----

# Read in impervious surface layer  
land <- rast("../Data/Plants/Plot Level/Berlin_impervious.tif")
plot(land)
land$Berlin_impervious

# Convert it to a df for plotting
land_df <- as.data.frame(land, xy = T)


# Read in garden sf
berlin_gardens <- sf::read_sf("../Data/Plants/Plot Level/Berlin_gardens.gpkg")


# Check crs match
st_crs(land)
st_crs(berlin_gardens)

# They dont, so convert gardens
berlin_gardens <- st_transform(berlin_gardens, crs = 3035)


# Plot

a <- ggplot() + 
  geom_raster(data = land_df, aes(x, y, fill = Berlin_impervious)) +  
  scale_fill_distiller(palette = "Oranges", direction = 1) + 
  geom_sf(data = berlin_gardens, 
          color = "black", pch = 16, size = 8.5) +   
  geom_sf(data = berlin_gardens, 
          color = "dodgerblue", pch = 16, size = 7) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        plot.title = element_text(size = 30, hjust = 0.5, face = "bold")) +
  labs(fill = "Impervious \nSurface Cover (%)") +
  ggtitle("Berlin")+
  ggspatial::annotation_scale()

a


#### Munich ----

# Read in Munich land raster
munich_land <- rast("../Data/Plants/Plot Level/Munich_impervious.tif")
munich_land_df <- as.data.frame(munich_land, xy = T)

# Now gardens
munich_gardens <- read_sf("../Data/Plants/Plot Level/Munich_gardens.gpkg")
munich_gardens <- st_transform(munich_gardens, crs = 3035)

# Plot
ggplot() + 
  geom_raster(data = munich_land_df, aes(x, y, fill = Munich_impervious)) +  
  scale_fill_distiller(palette = "Oranges", direction = 1) + 
  geom_sf(data = munich_gardens, 
          color = "black", pch = 16, size = 8.5) +   
  geom_sf(data = munich_gardens, 
          color = "dodgerblue", pch = 16, size = 7) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  labs(fill = "Impervious \nSurface Cover (%)")

# Remove Freising
munich_gardens_noF <- dplyr::filter(munich_gardens, gardenID != "KP")


b <- ggplot() + 
  geom_raster(data = munich_land_df, aes(x, y, fill = Munich_impervious)) +  
  scale_fill_distiller(palette = "Oranges", direction = 1) + 
  geom_sf(data = munich_gardens_noF, 
          color = "black", pch = 16, size = 8.5) +   
  geom_sf(data = munich_gardens_noF, 
          color = "dodgerblue", pch = 16, size = 7) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        plot.title = element_text(size = 30, hjust = 0.5, face = "bold")) +
  labs(fill = "Impervious \nSurface Cover (%)") +
  ggtitle("Munich")+
  ggspatial::annotation_scale()

b


ggspatial::annotation_scale()



# Put them together
ggpubr::ggarrange(a, b,
                  ncol = 2, nrow = 1,
                  common.legend = T,
                  legend = "bottom")













rm(list = ls())





