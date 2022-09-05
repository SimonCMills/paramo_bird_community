# Run and plot NMDS and run ANOSIM

# packages ----
library(dplyr); library(vegan); library(ggplot2)

# data ----
bird_data <- readRDS("data/paramo_bird_dataset.rds")
for_nmds <- bird_data %>%
    filter(!point == "CHP1") %>% # contains all-0 detections
    group_by(species, point) %>%
    ungroup %>%
    filter(!(species %in% c("Anas_andium", "Tringa_melanoleuca", "Spatula_discors"))) %>%
    select(species, point, Q, habitat_type) %>%
    unique %>%
    reshape2::dcast(., point ~ species, value.var="Q") 

# nmds ----
for_nmds2 <- for_nmds %>%
    select(-point)

nmds_out <- metaMDS(for_nmds2, k=4, trymax=200, distance="jaccard")

nmds_scores <- scores(nmds_out)$sites %>%
    as_tibble %>%
    bind_cols(point = for_nmds$point, .)  %>%
    left_join(., select(bird_data, point, habitat_type, elev_ALOS) %>% unique)

# plot nmds ----
cols <- c("1" = "#0072B2", "-1" = "#ee6a50", "Forest" = "#008763")
cols2 <- scales::alpha(cols, .5)

png("nmds_plot.png", height=100, width=170, units="mm", res=300)
par(mfrow=c(1,2))
# dims 1&2
plot(nmds_out, type="n", choices=c(1,2))
points(nmds_out, choices=c(1,2), display="sites", pch=21, bg=cols2[factor(nmds_scores$habitat_type)])
ordiellipse(nmds_out, choices=c(1,2), 
            groups=nmds_scores$habitat_type, col=cols[factor(nmds_scores$habitat_type)], conf=.95)

# dims1&3
plot(nmds_out, type="n", choices=c(3,2))
points(nmds_out, choices=c(3,2), display="sites", pch=21, bg=cols2[factor(nmds_scores$habitat_type)])
ordiellipse(nmds_out, choices=c(3,2), 
            groups=nmds_scores$habitat_type, col=cols[factor(nmds_scores$habitat_type)], conf=.95)
dev.off()

# anosim ----
habitat_info <- left_join(for_nmds, bird_data %>% select(point, habitat_type) %>% unique)
anosim_out <- anosim(for_nmds2, habitat_info$habitat_type, distance="bray")
anosim_out # p = 0.035; R = 0.14
saveRDS(anosim_out, "outputs/anosim_fit.rds")
