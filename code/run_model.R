# fit models
library(dplyr); library(flocker)

# datasets ----
bird_data <- readRDS("data/paramo_bird_dataset.rds")

# vegetation data
veg_data <- bird_data %>%
    select(point, elev_ALOS, habitat_type, site, cluster, 
           abu_espeletia, abu_shrub) %>%
    unique %>%
    mutate(elev_sc = scale(elev_ALOS))

# data summary ----
bird_data %>%
    ungroup %>%
    select(habitat_type, point) %>%
    unique %>%
    group_by(habitat_type) %>%
    summarise(n())

bird_data %>%
    ungroup %>%
    select(species) %>%
    unique %>%
    nrow


# format data for occupancy model ----
det <- select(bird_data, v1:v4) %>%
    mutate_all(.funs = function(x)ifelse(x==-99, NA, x)) %>% 
    as.matrix

cov <- select(bird_data, species, phylo, 
              abu_espeletia, abu_shrub, 
              elev_ALOS, relev, site, cluster, habitat_type, point) %>%
    mutate(sp_cl = interaction(species, cluster), 
           abu_espeletia_sc = scale(abu_espeletia), 
           abu_shrub_sc = scale(abu_shrub),
           sp_site = interaction(species, site), 
           elev_sc = scale(elev_ALOS), 
           habitat_sc = ifelse(habitat_type == "A", 1, -1))

time <- select(bird_data, hps1:hps4) %>%
    as.matrix
obsvr <- select(bird_data, obs1:obs4) %>%
    as.matrix

obsvr_species <- bird_data %>%
    mutate(obs1 = interaction(obs1, species), 
           obs2 = interaction(obs2, species), 
           obs3 = interaction(obs3, species), 
           obs4 = interaction(obs4, species)) %>%
    select(obs1:obs4) %>%
    as.matrix
vis_cov <- list(time = time, obsvr = obsvr, obsvr_species = obsvr_species)

# format for flocker
fd <- flocker::make_flocker_data(det, cov, vis_cov)

# functions
f_occ <- ~ 1 + elev_sc + (1|sp_site) + site + abu_espeletia_sc + 
    abu_shrub_sc + (1 + abu_espeletia_sc + abu_shrub_sc + elev_sc|species)

f_det <- ~ 1 + abu_espeletia_sc + abu_shrub_sc + time + (1 + time |species) + 
    obsvr

# fit occupancy ----
flocker_fit <- flock(f_occ, f_det, fd, rep_constant = F, 
                     backend = "cmdstanr", cores = 4, chains=4, 
                     file = "outputs/flocker_fit2.rds")

# fit vegetation ----
fit_esp <- brms::brm(bf(abu_espeletia ~ 1 + elev_sc + habitat_sc + site, 
                        zi ~ elev_sc), 
                     veg_data, adapt_delta = .99,
                     backend = "cmdstanr", 
                     family="zero_inflated_negbinomial",
                     chains = 4, cores = 4, 
                     file = "outputs/espeletia_fit.rds",
                     prior = c(set_prior("normal(0, 3)", "b")))

fit_shr <- brms::brm(bf(abu_shrub ~ 1 + elev_sc + habitat_sc + site, 
                        zi ~ elev_sc), 
                     veg_data, adapt_delta = .99,
                     backend = "cmdstanr", 
                     family="zero_inflated_negbinomial",
                     chains = 4, cores = 4, 
                     file = "outputs/shrub_fit.rds",
                     prior = c(set_prior("normal(0, 3)", "b")))
