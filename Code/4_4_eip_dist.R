# plot prevalences over time
# fit gam model to check for symmetruc logistic growth - not symmetric
# fit non-linear logistic bayesian model to prevalence data
# plot posterior prevalence over time (eip cdfs)
# plot summary of bayesian model (max prev, eip50, eip90-eip10)
# plot against survival, and plot area under the curve of their product (number of infectious days)
# plot areas under those curves as scatter plot ± iqr

# author: ivan casas

##################################################################################################################################
# Environment
##################################################################################################################################
source("Code/compile_data_survival.r") # sets the right wd and calls view() on each dataset
source("Code/functions.r")

# for interpreatbility, let's make ki270 the reference level...
surv$exposed <- as.factor(surv$exposed)
surv$temp_mean <- as.factor(surv$temp_mean)
surv$temp_range <- as.factor(surv$temp_range)
surv$species <- as.factor(surv$species)

surv$exposed <- relevel(surv$exposed, ref = "Control")
surv$temp_mean <- relevel(surv$temp_mean, ref = "27")
surv$temp_range <- relevel(surv$temp_range, ref = "0")
surv$species <- relevel(surv$species, ref = "An. gambiae")


# load packages
packages <- c(
    "ggplot2",
    "tidyr",
    "tidybayes",
    "emmeans",
    "gratia",
    "cowplot",
    "mgcv",
    "nlme",
    "flexsurv",
    "dplyr",
    "tidyverse",
    "betareg",
    "brms",
    "bayestestR",
    "GGally",
    "ggh4x",
    "lme4",
    "car"
)
for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}

##################################################################################################################################
# EIP distribution by treatment (from sporozoite prevalence data)
##################################################################################################################################
sp_prev <- surv |>
    filter(exposed == "Exposed") |>
    group_by(pot_unique, age, species, treatment, temp_mean, temp_range, replicate) |>
    summarise(
        mean_prevalence = mean(label_pf == "with Plasmodium", na.rm = TRUE), # proportion of positives among those qpcr'ed
        n_qpcr = sum(!is.na(label_pf)),
        mean_count = mean(pf_copies[pf_copies != 0], na.rm = TRUE) # mean dna count among all positives
    ) |>
    ungroup()
str(sp_prev)
view(sp_prev)

# quick preview
ggplot(sp_prev, aes(x = age, y = mean_prevalence, color = treatment, shape = replicate)) +
    geom_point() +
    scale_shape_manual(values = 1:nlevels(sp_prev$replicate)) +
    facet_grid(treatment ~ species) +
    # geom_smooth(
    #     method = "loess",
    #     se = FALSE
    # )
    geom_smooth(
        method = "glm",
        method.args = list(
            family = quasibinomial(link = "logit")
        ),
        se = FALSE
    )

# stats
##################################################################################################################################
surv_eip <- filter(surv, surv$exposed == "Exposed", !is.na(surv$label_pf)) # only with qPCR'ed exposed mosquitoes

str(surv_eip)
surv_eip$species <- as.factor(surv_eip$species)
surv_eip$temp_mean <- as.factor(surv_eip$temp_mean)
surv_eip$temp_range <- as.factor(surv_eip$temp_range)
surv_eip$treatment <- as.factor(surv_eip$treatment)
surv_eip$pot <- as.factor(surv_eip$pot)

surv_eip$age <- as.integer(surv_eip$age)

contrasts(surv_eip$species) <- contr.sum(nlevels(surv_eip$species))
contrasts(surv_eip$temp_mean) <- contr.sum(nlevels(surv_eip$temp_mean))
contrasts(surv_eip$temp_range) <- contr.sum(nlevels(surv_eip$temp_range))
contrasts(surv_eip$pot) <- contr.sum(nlevels(surv_eip$pot))


# First GAM to test for linearity of predictors with link function
# (i.e. interrogate symmetry of logistic growth for each group)

table(surv_eip$species, surv_eip$temp_mean, surv_eip$temp_range) # sample sizes are alright
surv_eip$grp <- interaction(surv_eip$species, surv_eip$temp_mean, surv_eip$temp_range)

m_gam <- gam(
    label_pf ~ s(age, by = grp) +
        species * temp_mean * temp_range +
        s(pot, bs = "re") + s(replicate, bs = "re"),
    family = binomial(link = "logit"),
    data = surv_eip,
    method = "REML" # (wood, 2011), also author of mgcv pkg
)
gam.check(m_gam) # tried with k = 20 and no difference. keeping low (default = 10) for speed
summary(m_gam)
# look at smooth effects (edf)
# edf is representative of the degree of complexity to expect (>1: nonlinear)


#  most groups show nonlinear effect of age on eip cdf  (linear logistic not supported):
#             edf     p
# s(age):G270   3    **
# s(age):C270   4   ***
# s(age):G210   2   ***
# s(age):C210   1     *
# s(age):G276   2    **
# s(age):C276   2
# s(age):G216   3    **
# s(age):C216   2     *
# s(pot).       0           pot isn't explaining much, can be dropped
# s(replicate)  6           replicate is an important random effect, must remain

# plot results
gam_eip <- draw(m_gam, select = 1:8, ncol = 4)
gam_eip

ggsave(plot = gam_eip, filename = "Figures/GAM_eip.png", width = 18, height = 12)
ggsave(plot = gam_eip, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/GAM_eip.png", width = 16, height = 12)

# conclusion:
# most cases show nonlinear growth
# no evidence for decay (only with sparse tail data)
# general evidence for slowed down growth on tail

# we'll include the random effect only on 1 parameter of the final model...
# let's decide where:

# random effect on max prev
bform_upper <- bf(
    label_pf ~ inv_logit(logitupper) * inv_logit(slope * (log(age) - log(ed50))),
    logitupper ~ 0 + grp + (1 | replicate),
    slope ~ 0 + grp,
    ed50 ~ 0 + grp,
    nl = TRUE
)

# random effect on eip50
bform_ed50 <- bf(
    label_pf ~ inv_logit(logitupper) * inv_logit(slope * (log(age) - log(ed50))),
    logitupper ~ 0 + grp,
    slope ~ 0 + grp,
    ed50 ~ 0 + grp + (1 | replicate),
    nl = TRUE
)

# random effect on grotwth rate
bform_slope <- bf(
    label_pf ~ inv_logit(logitupper) * inv_logit(slope * (log(age) - log(ed50))),
    logitupper ~ 0 + grp,
    slope ~ 0 + grp + (1 | replicate),
    ed50 ~ 0 + grp,
    nl = TRUE
) # doesnt converge well

# priors
priors <- c(
    prior(normal(2, 1.5), nlpar = "logitupper"), # logit(0.88) ≈ 2, expect asymptote near 70-95%
    prior(normal(2, 1.5), nlpar = "slope", lb = 0), # lb: lower bound
    prior(normal(20, 10), nlpar = "ed50", lb = 0) # eip50 prior
)

# fit them - this steps takes ages, only do it once...
m_upper <- brm(bform_upper,
    data = surv_eip, family = bernoulli(link = "identity"),
    prior = priors, chains = 4, cores = 4, iter = 4000,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)
m_ed50 <- update(m_upper, bform_ed50) # reuses compiled Stan code where possible
m_slope <- update(m_upper, bform_slope)

# not nested, and aic doesnt quite apply here cause random df <1. Compare with loo:
# also takes a minute
m_upper <- add_criterion(m_upper, "loo")
m_ed50 <- add_criterion(m_ed50, "loo")
m_slope <- add_criterion(m_slope, "loo") # convergence problems

comp <- loo_compare(m_upper, m_ed50)
print(comp, digits = 2)
z <- comp[2, "elpd_diff"] / comp[2, "se_diff"]
2 * pnorm(-abs(z)) # 0.07 - m_ed50 is better but it doesn't really matter that much really

# so, we specify our nonlinear logistic model
best_model <- m_ed50

# check model fit
summary(best_model)
# RHAt should be 1 (defo not >1) - good
# BulkESS and Tail ESS in the thousands - good

# make predictions and report effects as predicted differences
draws <- as_draws_df(best_model)
groups <- names(draws)[str_detect(names(draws), "^b_ed50_grp")] %>% str_remove("^b_ed50_grp")

get_param <- function(param, grp) {
    x <- draws[[paste0("b_", param, "_grp", grp)]]
    if (param == "logitupper") plogis(x) else x
}

# pull all 3 parameters for all 8 groups into one long df, one row per draw x group
post <- map_dfr(groups, function(g) {
    slope <- get_param("slope", g)
    ed50 <- get_param("ed50", g)
    upper <- get_param("logitupper", g)
    tibble(
        grp = g,
        draw = seq_along(slope),
        asymptote = upper,
        ed50 = ed50,
        eip90_eip10 = ed50 * (9^(1 / slope) - 9^(-1 / slope)), # usual dispersion metric
        mass_80 = upper * 0.8, # so eip90-10 can be scaled to the amount of events happening in that time interval
        iqr = ed50 * (3^(1 / slope) - 3^(-1 / slope)), # iqr
        mass_50 = upper * 0.5 # so iqr can be scaled to the amount of events happening in that time interval
    )
})

# which comparisons are relevant
contrasts <- list(
    "species at 21±0" = c("An.gambiae.21.0", "An.coluzzii.21.0"),
    "species at 21±6" = c("An.gambiae.21.6", "An.coluzzii.21.6"),
    "species at 27±0" = c("An.gambiae.27.0", "An.coluzzii.27.0"),
    "species at 27±6" = c("An.gambiae.27.6", "An.coluzzii.27.6"),
    "Tmean for AC ±0" = c("An.coluzzii.21.0", "An.coluzzii.27.0"),
    "Tmean for AC ±6" = c("An.coluzzii.21.6", "An.coluzzii.27.6"),
    "Tmean for AG ±0" = c("An.gambiae.21.0", "An.gambiae.27.0"),
    "Tmean for AG ±6" = c("An.gambiae.21.6", "An.gambiae.27.6"),
    "Trange for AC 21" = c("An.coluzzii.21.0", "An.coluzzii.21.6"),
    "Trange for AC 27" = c("An.coluzzii.27.0", "An.coluzzii.27.6"),
    "Trange for AG 21" = c("An.gambiae.21.0", "An.gambiae.21.6"),
    "Trange for AG 27" = c("An.gambiae.27.0", "An.gambiae.27.6")
)

# helpre to build table
summarise_diff <- function(x) {
    tibble(
        Median  = median(x),
        CI_low  = unname(quantile(x, .025)),
        CI_high = unname(quantile(x, .975)),
        pd      = max(mean(x > 0), mean(x < 0))
    )
}

table2 <- map_dfr(names(contrasts), function(cname) {
    g1 <- contrasts[[cname]][1]
    g2 <- contrasts[[cname]][2]
    d1 <- post %>% filter(grp == g1)
    d2 <- post %>% filter(grp == g2)
    map_dfr(c("asymptote", "ed50", "eip90_eip10", "mass_80", "iqr", "mass_50"), function(p) {
        diff <- d1[[p]] - d2[[p]]
        bind_cols(Comparison = cname, Parameter = p, summarise_diff(diff))
    })
}) %>%
    mutate(across(c(Median, CI_low, CI_high, pd), ~ round(., 3))) %>%
    arrange(Parameter, Comparison)

print(table2, n = Inf)

# aaaaand plot things
################################################################################################################
set.seed(1984)

selected_colors <- c("#002fff", "#80b5ff", "#ff0000", "#ff87eb")

# prediction grid
newdata <- expand.grid(
    age = seq(1, 30, length.out = 30), # we are gonna be doing log(age) later
    species = levels(surv_eip$species),
    temp_mean = levels(surv_eip$temp_mean),
    temp_range = levels(surv_eip$temp_range)
) |> mutate(
    treatment = paste0(temp_mean, "±", temp_range, "°C"),
    grp = interaction(species, temp_mean, temp_range, drop = TRUE),
    row_id = row_number()
)

# EIP predictions from bayesian model
preds <- posterior_epred(best_model, newdata = newdata, re_formula = NA)
B <- nrow(preds)

newdata$predicted <- apply(preds, 2, median)
newdata$lower <- apply(preds, 2, quantile, probs = 0.25)
newdata$upper <- apply(preds, 2, quantile, probs = 0.75)

# survival predictions from flexsurv model
surv_exp <- filter(surv, exposed == "Exposed") # only with exposed mosquitoes
best_dist <- "genf" # from paraemtric survival modelling script
m1 <- flexsurvreg(event ~ species * temp_mean * temp_range, data = surv_exp, dist = best_dist) # from paraemtric survival modelling script

nd1 <- expand.grid(
    species = levels(surv_eip$species),
    temp_range = levels(surv_eip$temp_range),
    temp_mean = levels(surv_eip$temp_mean)
)

pred1 <- summary(m1,
    newdata = nd1, type = "survival", ci = TRUE, tidy = TRUE,
    t = unique(newdata$age)
) %>%
    rename(age = time) %>% # will be useful later
    mutate(
        treatment = paste0(temp_mean, "±", temp_range, "°C"),
        S_t = est,
        est_logit = qlogis(est),
        se_logit = (qlogis(ucl) - qlogis(lcl)) / (2 * 1.96)
    )

# intersection of both predictions
area_df <- newdata %>%
    inner_join(pred1, by = c("age", "species", "treatment", "temp_mean", "temp_range")) %>%
    mutate(area_product = predicted * S_t)

nrow(newdata)
nrow(pred1)
nrow(area_df) # sanity check: should all match

# EIP posterior curves
################################################################################################################

eip_dists <- ggplot(newdata, aes(x = age, y = predicted, color = treatment, fill = treatment)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
    facet_wrap2(species ~ ., scales = "free") +
    scale_color_manual(values = selected_colors) +
    scale_fill_manual(values = selected_colors) +
    ylim(0, 0.8) +
    labs(
        x = "Days post exposure (dpe)",
        y = "Posterior prevalence (%) & IQR",
        colour = "Temperature",
        fill = "Temperature"
    ) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 35),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(colour = "black", linewidth = 0.6),
        plot.margin = margin(10, 10, 10, 10),
        strip.text = element_text(size = 35, face = "bold"),
        legend.position = "none",
    )
eip_dists
ggsave(plot = eip_dists, filename = "Figures/eip_dists.png", width = 18, height = 12)
ggsave(plot = eip_dists, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/eip_dists.png", width = 10, height = 7.5)

# extract legend for plotting
eip_dists_legend <- eip_dists + theme(
    legend.position = "right",
    legend.text = element_text(size = 15), # dont wanna change the theme in the plot
    legend.title = element_blank()
) +
    guides(
        color = guide_legend(ncol = 4, title.position = "top"),
        fill  = guide_legend(ncol = 4, title.position = "top")
    )

legend_only <- cowplot::get_legend(eip_dists_legend)
cowplot::save_plot(
    "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/colormap_2.png",
    legend_only,
    base_width = 5,
    base_height = 1
)

# posterior parameter values
################################################################################################################
# scale itself means nothing, so instead it's used to calculate eip90-eip10

lookup <- surv_eip %>% # each model names the same thing a different way sfruayibuoehiabvhk
    distinct(species, temp_mean, temp_range) %>%
    mutate(
        grp_clean = str_remove_all(interaction(species, temp_mean, temp_range, drop = TRUE), " "),
        treatment = paste0(temp_mean, "±", temp_range)
    )

post_summary <- post %>% # from earlier, posteriors from the model
    pivot_longer(c(asymptote, ed50, eip90_eip10), names_to = "Parameter", values_to = "value") %>%
    group_by(grp, Parameter) %>%
    summarise(
        Median = median(value),
        CI_low = quantile(value, .25),
        CI_high = quantile(value, .75),
        .groups = "drop"
    ) %>%
    left_join(lookup, by = c("grp" = "grp_clean")) %>%
    mutate(Parameter = factor(Parameter,
        levels = c("asymptote", "ed50", "eip90_eip10"),
        labels = c("Competence (%)", "EIP50 (d)", "EIP90-EIP10 (d)")
    ))

post_params <- ggplot(post_summary, aes(x = treatment, y = Median, color = treatment)) +
    geom_pointrange(aes(ymin = CI_low, ymax = CI_high), linewidth = 1, size = 2) +
    facet_grid2(species ~ Parameter, scales = "free", independent = "all") + # facet_wrap doesnt allow side titles :(
    labs(x = "Treatment", y = "Posterior estimate (median & IQR)") +
    scale_color_manual(values = selected_colors) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        axis.text.x = element_text(size = 26, angle = 25, hjust = 1),
        axis.title = element_text(size = 35),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(colour = "black", linewidth = 0.6),
        plot.margin = margin(10, 10, 10, 10),
        strip.text = element_text(size = 30, face = "bold"),
        legend.position = "none",
    )
post_params
ggsave(plot = post_params, filename = "Figures/post_params.png", width = 16, height = 12)
ggsave(plot = post_params, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/post_params.png", width = 13, height = 9.75)

# both curves and the area under their product
################################################################################################################

km1 <- get_km_data(surv_exp, nd1)

# unify levels so that the legends behave
treatment_levels <- sort(unique(as.character(sp_prev$treatment)))
sp_prev$treatment <- factor(sp_prev$treatment, levels = treatment_levels)
newdata$treatment <- factor(newdata$treatment, levels = treatment_levels)
pred1$treatment <- factor(pred1$treatment, levels = treatment_levels)
km1$treatment <- factor(km1$treatment, levels = treatment_levels)
area_df$treatment <- factor(area_df$treatment, levels = treatment_levels)

surv_eips <- ggplot(sp_prev, aes(x = age, y = mean_prevalence, colour = factor(treatment), shape = replicate)) +
    geom_point(size = 3) + # sp prevalence points empirical
    geom_line( # parametric prevalence (eip CDF)
        data = newdata, aes(x = age, y = predicted, color = factor(treatment)),
        inherit.aes = FALSE, linewidth = 1
    ) +
    scale_shape_manual(values = 1:nlevels(sp_prev$replicate)) +
    geom_line( # parametric survival
        data = pred1,
        aes(x = age, y = S_t, colour = factor(treatment)),
        linewidth = 1, inherit.aes = FALSE
    ) +
    geom_step( # empirical KM: step function
        data = km1,
        aes(x = time, y = est, colour = factor(treatment)),
        linewidth = 1, inherit.aes = FALSE
    ) +
    geom_ribbon(
        data = area_df,
        aes(x = age, ymin = 0, ymax = area_product, fill = factor(treatment)),
        inherit.aes = FALSE, alpha = 0.25
    ) +
    # facet_wrap2(species ~ treatment,
    #     scales = "free", ncol = 4,
    #     axes = "margins", remove_labels = "all"
    # ) +
    facet_grid2(species ~ treatment,
        scales = "free",
        axes = "margins", remove_labels = "all",
    ) +
    scale_color_manual(values = selected_colors) +
    scale_fill_manual(values = selected_colors) +
    scale_y_continuous(
        name = "Proportion alive at time t, S(t)",
        sec.axis = sec_axis(~., name = "Proportion infectious at time t, Prev(t)")
    ) +
    labs(
        x = "Age",
        colour = "Temperature",
        fill = "Temperature",
        shape = "Replicate"
    ) +
    guides(
        shape = guide_legend(order = 1, ncol = 5, title.position = "top"),
        color = guide_legend(order = 2, ncol = 2, title.position = "top"),
        fill  = guide_legend(order = 2, ncol = 2, title.position = "top")
    ) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        strip.text = element_text(size = 30, face = "bold"),
        strip.placement = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        plot.caption = element_text(size = 27, hjust = 0.5),
        axis.title = element_text(size = 30),
        plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 45, hjust = 0.5),
        legend.title = element_text(size = 30, hjust = 0.5),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.text = element_text(size = 24),
        legend.key.size = unit(1.5, "cm"),
    )
surv_eips
ggsave(plot = surv_eips, filename = "Figures/surv_eip.png", width = 18, height = 12)
ggsave(plot = surv_eips, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/surv_eip.png", width = 20, height = 14)


# wee sensitivity check for the extrapolation of KMs
################################################################################################################

last_death <- surv_exp %>%
    filter(dead == 1) %>%
    group_by(species, temp_mean, temp_range) %>%
    summarise(max_death_day = max(age))

pred1 %>%
    left_join(last_death, by = c("species", "temp_mean", "temp_range")) %>%
    group_by(species, temp_mean, temp_range) %>%
    summarise(pct_mass_extrapolated = sum(S_t[age > max_death_day]))

area_df_truncated <- area_df %>%
    left_join(last_death, by = c("species", "temp_mean", "temp_range")) %>%
    filter(age <= max_death_day) %>%
    group_by(species, treatment) %>%
    summarise(area_truncated = sum(area_product))

area_df <- area_df %>%
    group_by(species, treatment) %>%
    summarise(area_df_truncated_full = sum(area_product)) |>
    left_join(area_df_truncated, by = c("species", "treatment")) |>
    mutate(
        diff = abs(area_df_truncated_full - area_truncated), # amount of area from extrapolated survival
        proportion = diff / area_df_truncated_full
    ) # proportion of total area

mean(area_df$proportion)


# infectious days
################################################################################################################

pred1_aligned <- newdata %>%
    select(row_id, age, species, treatment, temp_mean, temp_range) %>%
    inner_join(pred1, by = c("age", "species", "treatment", "temp_mean", "temp_range"))

stopifnot(nrow(pred1_aligned) == nrow(newdata))

# simulate S_t draws on the logit scale
S_t_draws <- sapply(seq_len(nrow(pred1_aligned)), function(i) {
    plogis(rnorm(B, mean = pred1_aligned$est_logit[i], sd = pred1_aligned$se_logit[i]))
})

# elementwise product
area_draws <- preds[, pred1_aligned$row_id] * S_t_draws

# area under the product
groups_tbl <- newdata %>%
    distinct(species, treatment, row_id) %>%
    group_by(species, treatment) %>%
    summarise(row_ids = list(row_id), .groups = "drop")

auc_summary <- pmap_dfr(groups_tbl, function(species, treatment, row_ids) {
    idx <- match(row_ids, newdata$row_id)

    auc_draws <- rowSums(area_draws[, idx, drop = FALSE]) # discrete sum, per posterior draw

    tibble(
        species = species,
        treatment = treatment,
        AUC_median = median(auc_draws),
        AUC_low = quantile(auc_draws, 0.25), # 0.025
        AUC_high = quantile(auc_draws, 0.75) # 0.975
    )
})

print(auc_summary, n = Inf)


auc_plot <- ggplot(auc_summary, aes(x = treatment, y = AUC_median, color = treatment)) +
    geom_point(size = 8) +
    geom_errorbar(aes(ymin = AUC_low, ymax = AUC_high), width = 0.15, linewidth = 0.8) +
    facet_wrap(~species, scales = "free_x") +
    scale_color_manual(values = selected_colors) +
    labs(x = NULL, y = "Area under predicted infectious-prevalence curve\n(IQR)") +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 35),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(colour = "black", linewidth = 0.6),
        plot.margin = margin(10, 10, 10, 10),
        strip.text = element_text(size = 35, face = "bold"),
        legend.position = "none",
    )

auc_plot
ggsave(plot = auc_plot, filename = "Figures/inf_days.png", width = 18, height = 12)
ggsave(plot = auc_plot, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/inf_days.png", width = 16, height = 12)
