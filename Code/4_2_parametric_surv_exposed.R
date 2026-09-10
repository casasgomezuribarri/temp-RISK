# survival analysis ch4
# fit cox model - invalid
# fit parametric model - genF
# plot parametric fit + km

# author: ivan casas


##################################################################################################################################
# Environment (loading packages and data)
##################################################################################################################################

source("Code/compile_data_survival.r")
source("Code/functions.r")

# load packages
packages <- c(
    "knitr",
    "dplyr",
    "survival",
    "emmeans",
    "ggplot2",
    "tibble",
    "devtools",
    "readr",
    "lubridate",
    "DT",
    "ggsurvfit",
    "gtsummary",
    "tidycmprsk",
    "RColorBrewer",
    "survminer",
    "coxme",
    "rms",
    "gridExtra",
    "flexsurv",
    "muhaz",
    "data.table",
    "readxl",
    "gridExtra"
)
for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}

surv$exposed <- as.factor(surv$exposed)
surv$species <- as.factor(surv$species)
surv$temp_mean <- as.factor(surv$temp_mean)
surv$temp_range <- as.factor(surv$temp_range)
surv$treatment <- as.factor(surv$treatment)
surv$age <- surv$age
str(surv)
surv_exp <- filter(surv, exposed == "Exposed") # only with exposed mosquitoes

##################################################################################################################################
# cox model
##################################################################################################################################

# simplest model (no strat, no interactions)
full <- coxph(event ~ species + temp_range + temp_mean + cluster(replicate), data = surv_exp) # always inspect with/without replicate
summary(full)
survival::cox.zph(full) # invalid model

# Schoenfeld residuals against time - in a good model these lines would be straight
plot(cox.zph(full), var = "species")
plot(cox.zph(full), var = "temp_range")
plot(cox.zph(full), var = "temp_mean")

##################################################################################################################################
# parametric survival: only exposed, effect of temp and species
##################################################################################################################################
# flexsurvreg does not accept random effects

# choose a distribution
par_fit <- compare_parametric_fits(
    data = surv_exp,
    time_var = "age",
    event_var = "dead",
    plot_title = paste0("Parametric fits")
)
# save plot
par_fit$plot
png(
    file = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/parafits_exposed.png",
    width = 800, height = 800
)
par_fit$plot
dev.off()

par_fit$comparison

best_dist <- "genf"

m1 <- flexsurvreg(event ~ species * temp_mean * temp_range, data = surv_exp, dist = best_dist)
report <- tidy(m1)

print(report)
view(report)

# p computation - type iii analysis framework. This is not supported in car pkg, so we do it manually:
##################################################################################################################################

# we'll need to reencode variables to sum-to-0
local_contrasts <- list(
    species = contr.sum(2),
    temp_mean = contr.sum(2),
    temp_range = contr.sum(2)
)

# using those contrasts, build the model matrix manually
X_matrix <- model.matrix(~ species * temp_mean * temp_range,
    data = surv_exp,
    contrasts.arg = local_contrasts
)[, -1]

# gotta remove colons from col names
colnames(X_matrix) <- gsub(":", "_", colnames(X_matrix))

# bind to surv data
surv_data_sumto0 <- cbind(surv_exp, as.data.frame(X_matrix))


# full model formula (model.matrix warped names a bit but thats useful)
full_formula <- event ~ species1 + temp_mean1 + temp_range1 + # mains
    species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + # 2ways
    species1_temp_mean1_temp_range1 # 3way

# this should be equivalent to the og model
best_model_sumto0 <- flexsurvreg(full_formula, data = surv_data_sumto0, dist = best_dist) # used to be full_surv

# sanity check that they are indeed the same model
logLik(m1) # -8108.464 (df=11)
logLik(best_model_sumto0) # -8108.464 (df=11)
AIC(m1) # 16238.93
AIC(best_model_sumto0) # 16238.93

# right. now fit every single relevant nested model (full model without 1 term)
no_sp <- flexsurvreg(event ~ temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_tm <- flexsurvreg(event ~ species1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_tr <- flexsurvreg(event ~ species1 + temp_mean1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_sp_tm <- flexsurvreg(event ~ species1 + temp_mean1 + temp_range1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_sp_tr <- flexsurvreg(event ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_tm_tr <- flexsurvreg(event ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_sp_tm_tr <- flexsurvreg(event ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)

# custom function to show them all together
p_table_surv <- rbind(
    calc_lrt(best_model_sumto0, no_sp, "species"),
    calc_lrt(best_model_sumto0, no_tm, "temp_mean"),
    calc_lrt(best_model_sumto0, no_tr, "temp_range"),
    calc_lrt(best_model_sumto0, no_sp_tm, "species:temp_mean"),
    calc_lrt(best_model_sumto0, no_sp_tr, "species:temp_range"),
    calc_lrt(best_model_sumto0, no_tm_tr, "temp_mean:temp_range"),
    calc_lrt(best_model_sumto0, no_sp_tm_tr, "species:temp_mean:temp_range")
)

print(p_table_surv)
view(report) # coefs

# visualise predictions
##################################################################################################################################

nd1 <- bind_rows(
    expand.grid(
        species = c("An. gambiae", "An. coluzzii"),
        temp_range = as.factor(c(0, 6)),
        temp_mean = as.factor(c(21, 27))
    )
)

# add counts for plot annotations
nd1$count <- numeric(nrow(nd1))
for (i in 1:nrow(nd1)) {
    # pot info
    sp <- nd1$species[i]
    tm <- nd1$temp_mean[i]
    tr <- nd1$temp_range[i]

    # grab pot data from surv
    pot <- subset(surv_exp, species == sp & temp_mean == tm & temp_range == tr)
    # count rows
    nd1$count[i] <- nrow(pot)
}

# sanity check - these should be the same
sum(nd1$count)
nrow(surv_exp)

# make predictions on the new synthetic datasets
pred1 <- summary(m1, newdata = nd1, type = "survival", ci = TRUE, tidy = TRUE)

pred1 <- pred1 %>%
    mutate(
        treatment = paste0(temp_mean, "±", temp_range, "°C")
    )


# kms to overlap with predictions (custom function for computing kms)
km1 <- get_km_data(surv_exp, nd1)

# calculate median survival per facet and group
median_survival <- km1 %>%
    group_by(species, treatment) %>%
    summarise(
        median_t = min(time[est <= 0.5], na.rm = TRUE),
        .groups = "drop"
    )

# add y positions for plotting
median_annotations <- median_survival %>%
    group_by(species) %>% # in each facet there's 4 curves
    arrange(species, treatment) %>% # give them idx
    mutate(
        y_pos = case_when(
            row_number() == 1 ~ 0.90,
            row_number() == 2 ~ 0.84,
            row_number() == 3 ~ 0.78,
            TRUE ~ 0.72 # last one (TRUE is 'all the rest')
        )
    ) %>%
    ungroup() %>%
    mutate(label = paste0(median_t))

# probs good to add title to annotation
median_title <- median_annotations %>%
    distinct(species, treatment) %>%
    mutate(
        y_pos = 0.97,
        label = "Median (d)"
    )

selected_colors <- c("#002fff", "#80b5ff", "#ff0000", "#ff87eb")

# plot preds + km + median annotations
parapreds1_km <- ggplot(pred1, aes(x = time, y = est, colour = factor(treatment))) +
    geom_line(linewidth = 1.5) + # parametric fitted lines
    geom_ribbon( # confidence intervals
        aes(ymin = lcl, ymax = ucl, fill = factor(treatment)),
        alpha = 0.15, colour = NA
    ) + # colour = NA to avoid border around ribbons
    geom_step( # empirical KM: step function, same colour mapping, no legend duplication
        data = km1,
        aes(x = time, y = est, colour = factor(treatment)),
        linewidth = 1.5, linetype = "dashed", inherit.aes = FALSE
    ) +
    geom_segment( # vertical lines from y=0.5 to y=0 at x=median_t
        data = median_survival,
        aes(
            x = median_t, xend = median_t,
            y = 0.5, yend = 0,
            colour = factor(treatment)
        ),
        linewidth = 0.8, linetype = "solid", inherit.aes = FALSE
    ) +
    geom_text( # 'median' titles (black font bold)
        data = median_title,
        aes(x = 20, y = y_pos, label = label),
        hjust = 0, vjust = 0,
        size = 8, # fontface = "bold",
        colour = "black",
        inherit.aes = FALSE
    ) +
    geom_text( # median values in respective colours
        data = median_annotations,
        aes(x = 25, y = y_pos, label = label, colour = factor(treatment)),
        hjust = 0, vjust = 0,
        size = 8, fontface = "bold",
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    facet_wrap(~species, scales = "free", ncol = 2) +
    scale_fill_manual(values = selected_colors) +
    scale_color_manual(values = selected_colors) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
        # title = "Empirical and Predicted Survival Curves",
        x = "Days post infection",
        y = "Survival probability",
        colour = "Temperature",
        fill = "Temperature",
        # caption = "Solid = parametric fit; dashed = Kaplan-Meier" # no caption in publication, this is for presentations
    ) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        plot.caption = element_text(size = 27, hjust = 0.5),
        strip.text = element_text(size = 35, face = "bold"),
        # strip.background = element_rect(fill = "#e8e7ff", colour = "#8c8cff", linewidth = 2),
        axis.title = element_text(size = 35),
        plot.margin = margin(10, 10, 10, 10),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(colour = "black", linewidth = 0.6),
        # plot.title = element_text(size = 45, hjust = 0.5), # not title for publication
        legend.position = "none", # no legend for publication (it's part of a panel)
        # legend.title = element_text(size = 30), # Font size for legend title
        # legend.text = element_text(size = 24), # Font size for legend text
        # legend.key.size = unit(1.5, "cm") # Size of legend keys
    )

parapreds1_km
ggsave(plot = parapreds1_km, "Figures/parametric_exposed_km_ann.png", width = 16, height = 12, units = "in", dpi = 150)
ggsave(plot = parapreds1_km, "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/parametric_km_only_exposed_ann.png", width = 10, height = 7.5, units = "in", dpi = 150)


# actual coefficients (output is in log scale)
coefs <- tidy(m1, conf.int = TRUE) %>%
    filter(!term %in% c("mu", "sigma", "Q", "P")) %>% # remove dist params
    mutate(
        TR = exp(estimate), # time ratio = exp(coef)
        TR_lo = exp(conf.low),
        TR_hi = exp(conf.high),
        term = factor(term, levels = rev(unique(term))) # make it a factor with levels in (rev) row order so plotting is in order (as opposed to alphabetically)
    )
view(coefs)

# plot parameter estimates
coefs_plot <- coefs %>%
    ggplot(aes(x = TR, y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = TR_lo, xmax = TR_hi), height = 0.2) +
    geom_point(size = 3, colour = "#1a6fa8") +
    scale_x_log10() +
    labs(x = "Time Ratio (log scale)", y = NULL) +
    theme_bw(base_size = 20)
coefs_plot

ggsave(plot = coefs_plot, "Figures/parametric_exposed_coefs.png", width = 15, height = 17, units = "in", dpi = 150)
ggsave(plot = coefs_plot, "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/parametric_exposed_coefs.png", width = 10, height = 12, units = "in", dpi = 150)
