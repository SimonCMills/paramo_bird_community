# Packages ----
library(dplyr)

# read in the full bird dataset & subset to paramo points only
bird_data <- readRDS("data/birds.RDS")

paramo_bird_data <- bird_data %>% 
    ungroup %>%
    filter(site %in% c("Chingaza", "Tutaza", "Belen") |
               point %in% c(paste0("PUP", 1:6), paste0("PUF", 4:6),
                            paste0("IGF", 10:12), paste0("IGP", 1:3)))

# Read in dataset and format ----
# vegetation data from Edicson
veg_data <- data.table::fread("data/Espeletia_paramo.csv") %>%
    mutate(point = paste0(toupper(site), point)) %>%
    select(-site) %>%
    mutate(point = case_when(point == "IGA1" ~ "IGF10", 
                             point == "IGA2" ~ "IGF11", 
                             point == "IGA3" ~ "IGF12", 
                             TRUE ~ point)) %>%
    mutate(Espeletia = ifelse(point == "CHA22", 0, Espeletia), 
           shrubs = ifelse(point == "CHA22", 22, shrubs)) %>% 
    rename(abu_espeletia = Espeletia, abu_shrub = shrubs)

species_lookup <- read.csv("../Colombia_rangeLims/data/initial_species_list.csv", as.is=T) %>%
    select(species = HBW, species_eltontraits = eltontraits, species_clements = eBird) %>%
    mutate_all(function(x) gsub(" ", "_", x))

# Vegetation dataset ----
# there are four points missing in this dataset: CHA10D-12D, and CHA22 (i.e. 
# all paramo points, not pasture). According to Edicson, CHA10D-12D never got 
# vegetation done, and CHA22 is 0ESP 22SHR. 
pbd2 <- paramo_bird_data %>% 
    left_join(., species_lookup) %>%
    left_join(., veg_data) %>%
    mutate(habitat_type = substr(point, 3, 3), 
           habitat_type = case_when(point %in% paste0("BEA", 10:12) ~ "P", 
                                    point %in% paste0("TUA", 19:21) ~ "P", 
                                    point %in% paste0("IGF", 10:12) ~ "A", 
                                    point %in% paste0("PUF", 4:6) ~ "A",
                                    TRUE ~ habitat_type)) %>%
    filter(!habitat_type == "F") %>%
    filter(!grepl("PU", substr(point, 1, 2))) %>%
    mutate(site = case_when(site %in% c("Belen", "Tutaza") ~ "La Rusia", 
                            site == "IG" ~ "Iguaque", 
                            TRUE ~ site)) %>%
    filter(!(species %in% c("Anas_andium", "Tringa_melanoleuca", "Spatula_discors"))) %>%
    group_by(species) %>%
    filter(any(Q==1)) %>%
    ungroup %>%
    filter(!is.na(abu_espeletia))

df_veg <- pbd2 %>%
    select(point, abu_espeletia, abu_shrub, habitat_type, 
           elev_ALOS) %>%
    unique %>%
    mutate(habitat_sc = ifelse(habitat_type == "A", 1, -1))

saveRDS(df_veg, "data/paramo_vegetation_dataset.rds")
saveRDS(pbd2, "data/paramo_bird_dataset.rds")
