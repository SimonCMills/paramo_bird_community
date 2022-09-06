# Plot model results


# packages ----
library(brms); library(flocker); library(ggplot2); library(dplyr)

# plot theme
theme_eff_plot <- theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour = "black"))
cols <- c("A" = "#0072B2", "P" = "#ee6a50")

# predict vegetation ----
elev_range <- range(veg_data$elev_ALOS)
elev_scaling <- attr(veg_data$elev_sc, "scaled:scale")
elev_cent <- attr(veg_data$elev_sc, "scaled:center")

elev_scaling
pred_dat <- expand.grid(elev_ALOS = seq(elev_range[1], elev_range[2], len=100),
                        habitat_type = c("A", "P"),
                        site = "Chingaza",
                        cluster=NA, point=NA) %>%
    mutate(habitat_sc = ifelse(habitat_type == "P", -1, 1), 
           elev_sc = (elev_ALOS - elev_cent)/elev_scaling)

fits_shr <- fitted(fit_shr, re_formula = NULL, newdata = pred_dat,
                   summary = FALSE) %>%
    t()

fits_esp <- fitted(fit_esp, re_formula = NULL, newdata = pred_dat, 
                   summary = FALSE) %>%
    t()

# summarise fits ----
fits_shr_summ <- tibble(estimate = matrixStats::rowMeans2(fits_shr), 
                          lwr = matrixStats::rowQuantiles(fits_shr, probs = c(.05)), 
                          upr = matrixStats::rowQuantiles(fits_shr, probs = c(.95))) %>%
    bind_cols(pred_dat, .) %>%
    as_tibble

fits_esp_summ <- tibble(estimate = matrixStats::rowMeans2(fits_esp), 
                          lwr = matrixStats::rowQuantiles(fits_esp, probs = c(.05)), 
                          upr = matrixStats::rowQuantiles(fits_esp, probs = c(.95))) %>%
    bind_cols(pred_dat, .) %>%
    as_tibble


## plot ----
p1 <- ggplot(fits_esp_summ, aes(elev_ALOS, estimate)) + 
    geom_point(data=veg_data, aes(y=abu_espeletia, fill=habitat_type), pch=21, 
               alpha=.3, size=2) +
    geom_point(data=veg_data, aes(y=abu_espeletia), pch=21, col="black", size=2) +
    geom_line(aes(col=factor(habitat_type))) +
    geom_ribbon(aes(ymin=lwr, ymax = upr, fill=factor(habitat_type)), alpha=.2) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black"), 
          axis.text.x = element_blank()) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols) +
    guides(fill="none", colour="none") +
    labs(x = "", y = "Espeletia abundance")

p2 <- ggplot(fits_shr_summ, aes(elev_ALOS, estimate)) + 
    geom_point(data=veg_data, aes(y=abu_shrub, fill=habitat_type), pch=21, 
               alpha=.3, size=2) +
    geom_point(data=veg_data, aes(y=abu_shrub), pch=21, col="black", size=2) +
    geom_line(aes(col=factor(habitat_type))) +
    geom_ribbon(aes(ymin=lwr, ymax = upr, fill=factor(habitat_type)), alpha=.3) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black")) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols) +
    guides(fill="none", colour="none") +
    labs(x = "Elevation (m.a.s.l.)", y = "Shrub abundance")

p_both <- egg::ggarrange(p1, p2, ncol=1)
ggsave("figures/plot_vegetation.png", plot = p_both, width=100*1.2, 
       height = 150*1.2, dpi=400, units="mm")

# predictions ----
shr_center <- attributes(cov$abu_shrub_sc)$`scaled:center`
shr_scaling <- attributes(cov$abu_shrub_sc)$`scaled:scale`

esp_center <- attributes(cov$abu_espeletia_sc)$`scaled:center`
esp_scaling <- attributes(cov$abu_espeletia_sc)$`scaled:scale`

# predict occupancy 
n_draws <- 1000
draw_interval <- 4000 %/% n_draws
draw_ids <- seq(1, 4000, draw_interval)

shr_fit_vector <- fits_shr[,draw_ids] %>%
    as.vector

esp_fit_vector <- fits_shr[,draw_ids] %>%
    as.vector

pred_occupancy_df <- replicate(n_draws, pred_dat, FALSE) %>%
    bind_rows(., .id = "id_draw") %>%
    select(-habitat_type, -site) %>%
    mutate(draw = draw_ids[as.numeric(id_draw)], 
           abu_shrub = shr_fit_vector, 
           abu_espeletia = esp_fit_vector, 
           elev_sc = elev_sc,
           abu_shrub_sc = (abu_shrub - shr_center)/shr_scaling, 
           abu_espeletia_sc = (abu_espeletia - esp_center)/esp_scaling)

n_species <- length(unique(fd$data$species))
species <- dimnames(coefs_occ_sp)[[2]]
coefs_occ <- coef(flocker_fit, summary = F)
coefs_occ_sp <- coefs_occ$species

coef_occ_df <- tibble(species = rep(dimnames(coefs_occ_sp)[[2]], each = dim(coefs_occ_sp)[1]),
       draw = rep(1:dim(coefs_occ_sp)[1], dim(coefs_occ_sp)[2]),
       intercept = as.vector(coefs_occ_sp[,,"occ_Intercept"]), 
       occ_abu_espeletia_sc = as.vector(coefs_occ_sp[,,"occ_abu_espeletia_sc"]),
       occ_abu_shrub_sc = as.vector(coefs_occ_sp[,,"occ_abu_shrub_sc"]),
       occ_elev_sc = as.vector(coefs_occ_sp[,,"occ_elev_sc"])) %>%
    filter(draw %in% unique(pred_occupancy_df$draw))

occ_preds <- left_join(pred_occupancy_df, coef_occ_df) %>%
    mutate(p = boot::inv.logit(intercept + 
               occ_abu_espeletia_sc * abu_espeletia_sc + 
               occ_abu_shrub_sc * abu_shrub_sc + 
               occ_elev_sc * elev_sc))

# sense check occ ----
# occ_preds %>%
#     filter(species == "Spinus_spinescens") %>%
#     ggplot(aes(elev_sc, p, col = factor(habitat_sc), group=interaction(habitat_sc, draw))) + geom_line()
# 
# occ_preds %>%
#     group_by(elev_sc, habitat_sc, draw) %>%
#     summarise(SR = sum(p)) %>%
#     ggplot(aes(elev_sc, SR, col=factor(habitat_sc), group=interaction(habitat_sc, draw))) +
#     geom_line()

## SR plots ----
occ_summ <- occ_preds %>%
    group_by(elev_ALOS, habitat_sc, draw) %>%
    summarise(SR = sum(p)) %>%
    group_by(elev_ALOS, draw) %>%
    summarise(SR_P = SR[1], SR_A = SR[2], SR_diff = SR_P - SR_A) %>%
    group_by(elev_ALOS) %>%
    summarise(estimate_P = mean(SR_P), 
              lwr_P = quantile(SR_P, .05), 
              upr_P = quantile(SR_P, .95), 
              estimate_A = mean(SR_A), 
              lwr_A = quantile(SR_A, .05), 
              upr_A = quantile(SR_A, .95), 
              estimate_diff = mean(SR_diff), 
              lwr_diff = quantile(SR_diff, .05), 
              upr_diff = quantile(SR_diff, .95))

p1 <- ggplot(occ_summ, aes(elev_ALOS, estimate_A, ymin = lwr_A, ymax = upr_A)) + 
    geom_line(col = cols[1]) +
    geom_ribbon(alpha=.3, fill = cols[1]) +
    geom_line(aes(y = estimate_P), col=cols[2]) +
    geom_ribbon(aes(ymin = lwr_P, ymax = upr_P), alpha=.3, fill=cols[2]) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black"), 
          axis.text.x = element_blank()) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols) +
    guides(fill="none", colour="none") +
    labs(x = "", y = "Species richness")

p2 <- ggplot(occ_summ, aes(elev_ALOS, estimate_diff, ymin = lwr_diff, ymax = upr_diff)) + 
    geom_line(col = "black") +
    geom_ribbon(alpha=.1, fill = "black") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black")) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols) +
    geom_hline(yintercept = 0, lty = "longdash") +
    guides(fill="none", colour="none") +
    scale_y_continuous(breaks=seq(-20, 20, 2)) +
    labs(x = "Elevation (m.a.s.l.)", y = "Difference")

p_both <- egg::ggarrange(p1, p2, ncol=1, heights = c(1, .5))
ggsave("figures/plot_species_richness.png", plot = p_both, 
       width=100*1.2, height = 120*1.2, dpi=400, units="mm")

# coefficient plot ----
## fixed effects ----
feffs <- brms::fixef(flocker_fit, probs=c(.05, .95)) %>%
    as.data.frame %>%
    mutate(par = row.names(.)) %>%
    select(everything(), mid = Estimate, lwr = 3, upr=4) 

feffs_esp <- feffs %>%
    filter(par == "occ_abu_espeletia_sc")

feffs_shr <- feffs %>%
    filter(par == "occ_abu_shrub_sc")

# species-level effects ----
coefs_summ <- coef(flocker_fit)
occ_shr_summ <- coefs_summ$species[,,"occ_abu_shrub_sc"] %>%
    as.data.frame %>%
    mutate(species = gsub("_", " ", row.names(.))) %>%
    rename(estimate = Estimate, lwr = 3, upr = 4) %>%
    as_tibble %>%
    arrange(estimate) %>%
    mutate(species =factor(species, levels=species))

occ_esp_summ <- coefs_summ$species[,,"occ_abu_espeletia_sc"] %>%
    as.data.frame %>%
    mutate(species = gsub("_", " ", row.names(.))) %>%
    rename(estimate = Estimate, lwr = 3, upr = 4) %>%
    as_tibble %>%
    arrange(estimate) %>%
    mutate(species =factor(species, levels=levels(occ_shr_summ$species)))

#
xrange_shr <- c(min(occ_shr_summ$lwr), max(occ_shr_summ$upr))
xrange_esp <- c(min(occ_esp_summ$lwr), max(occ_esp_summ$upr))

p1 <- ggplot(occ_shr_summ, aes(estimate, species)) +
    geom_point() +
    geom_linerange(aes(xmin=lwr, xmax=upr)) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black"), 
          axis.text.x = element_blank(),
          axis.title = element_blank(), 
          plot.title = element_text(size = 11)) +
    scale_x_continuous(limits=xrange_shr, breaks = seq(-4, 4, 1)) +
    geom_vline(xintercept = 0, lty = "longdash") +
    guides(fill="none", colour="none") +
    labs(title = "(a) Effect of shrub abundance")

p2 <- ggplot(occ_esp_summ, aes(estimate, species)) +
    geom_point() +
    geom_linerange(aes(xmin=lwr, xmax=upr)) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black"), 
          axis.text.y = element_blank(), 
          axis.text.x = element_blank(),
          axis.title = element_blank(), 
          plot.title = element_text(size = 11)) +
    scale_x_continuous(limits=xrange_esp, breaks = seq(-4, 4, 1)) +
    geom_vline(xintercept = 0, lty = "longdash") +
    guides(fill="none", colour="none")  +
    labs(title = "(a) Effect of espeletia abundance")

# effects plots for shrubs
p3 <- ggplot(feffs_shr, aes(mid, "Average")) + geom_point() +
    geom_errorbarh(aes(xmin=lwr, xmax=upr), height=0) +
    geom_vline(xintercept = 0, lty="longdash") +
    theme_eff_plot +
    scale_x_continuous(limits=xrange_shr, breaks = seq(-4, 4, 1)) +
    labs(y = "", x = "Estimate (\u00B1 90% CI)")

p4 <- ggplot(feffs_esp, aes(mid, "Average")) + geom_point() +
    geom_errorbarh(aes(xmin=lwr, xmax=upr), height=0) +
    geom_vline(xintercept = 0, lty="longdash") +
    theme_eff_plot +
    theme(axis.text.y = element_blank()) +
    scale_x_continuous(limits=xrange_esp, breaks = seq(-4, 4, 1)) +
    labs(y = "", x = "Estimate (\u00B1 90% CI)")

p_both <- egg::ggarrange(p1, p2, p3, p4, ncol=2, heights=c(1, .05))
ggsave("figures/effect_plot_shrubs.png", p_both, units = "mm", 
       width = 170, height = 200)
