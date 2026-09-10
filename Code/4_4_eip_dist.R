# plot prevalences over time
# fit gam model to check for symmetruc logistic growth - not symmetric
# fit non-linear logistic bayesian model to prevalence data (priors from gamm)
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

surv$exposed <- as.factor(surv$exposed)
surv$temp_mean <- as.factor(surv$temp_mean)
surv$temp_range <- as.factor(surv$temp_range)
surv$species <- as.factor(surv$species)

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
    "bayesplot",
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
# prevalence from qpcr data
sp_prev <- surv |>
    filter(exposed == "Exposed") |> # only exposed samples
    group_by(pot_unique, dpi, age, species, treatment, temp_mean, temp_range, replicate) |>
    summarise(
        prevalence = mean(label_pf == "with Plasmodium", na.rm = TRUE), # proportion of positives among those qpcr'ed
        n = sum(!is.na(label_pf)),
        method = "qpcr"
    ) |>
    ungroup() |>
    filter(n > 0) # keep only rows that represent actual qpcr data
# str(sp_prev)
# table(surv$dead)
# view(sp_prev)

pot_meta <- surv |> # store matadata
    distinct(pot_unique, dpi, age, species, treatment, temp_mean, temp_range, replicate, exp_age)

# sporozoites observed through microscopy
spz_observed <- dissection |> # one row per mosquito
    group_by(pot_unique, dpi) |> # for each pot, rep, dpi
    filter(sum(!is.na(sporozoites)) > 0) |> # keep only samples from days where at least one gland was dissected
    summarise(
        n = sum(!is.na(sporozoites)), # total salivary glands dissected
        n_positive = sum(sporozoites, na.rm = TRUE), # total with spz
        prevalence = n_positive / n,
        method = "dissection"
    ) |>
    left_join(pot_meta, by = c("pot_unique", "dpi")) |> # add metadata stored in pot_meta
    select(
        pot_unique, dpi, age, species, treatment, temp_mean, temp_range,
        replicate, prevalence, n, method
    ) # same order as in sp_prev
# view(spz_observed)

combined_prev <- bind_rows(sp_prev, spz_observed) |>
    arrange(pot_unique, dpi) # one row is a compbination of pot, rep and dpi

# quick preview
ggplot(combined_prev, aes(x = age, y = prevalence, color = treatment, weight = n)) + # , shape = replicate
    geom_point(size = 3) +
    facet_grid(treatment ~ species) +
    xlim(c(0, 35)) +
    geom_smooth(
        method = "glm",
        method.args = list(
            family = quasibinomial(link = "logit")
        ),
        se = FALSE
    )

ggplot(combined_prev, aes(x = age, y = prevalence, color = treatment, shape = replicate, weight = n)) +
    geom_point(size = 5) +
    scale_shape_manual(values = 1:nlevels(sp_prev$replicate)) +
    facet_grid(treatment ~ species) +
    xlim(c(0, 35)) +
    geom_smooth(
        method = "glm",
        method.args = list(
            family = quasibinomial(link = "logit")
        ),
        se = FALSE
    )

ggplot(combined_prev, aes(x = age, y = prevalence, color = treatment, weight = n)) +
    geom_point(size = 5) +
    scale_shape_manual(values = 1:nlevels(sp_prev$replicate)) +
    facet_grid(replicate ~ species) +
    xlim(c(0, 35)) +
    geom_smooth(
        method = "glm",
        method.args = list(
            family = quasibinomial(link = "logit")
        ),
        se = FALSE
    )

# stats
##################################################################################################################################
# subset surv to calculate eip CDF
surv_eip <- surv |> # one row per mosquito
    filter(
        surv$exposed == "Exposed", # only exposed
        # !(surv$replicate %in% c("r00", "r01", "r02", "r03")), # considered removing replicates with qpcr but no dissection, but fit was worse
        !is.na(surv$label_pf) # only qpcr'ed
    ) |>
    select( # select the cols that are useulf for this analysis
        pot_unique, dpi, age, species, treatment, temp_mean, temp_range,
        replicate, pot, id, exp_age, label_pf
    ) |>
    mutate(method = "qpcr")

# add the dissecteds
spz_diss <- dissection |> # one row per mosquito
    filter(!is.na(sporozoites)) |> # include only mossies with dissected salivary glands
    left_join(pot_meta, by = c("pot_unique", "dpi", "replicate")) |> # add metadata stored in pot_meta so it it matches surv_eip. pit_unique, dpi, and replicate already present
    mutate(
        label_pf = case_when( # add label_pf based on whether sporozoites were observed for this mosquito
            sporozoites == TRUE ~ "with Plasmodium",
            sporozoites == FALSE ~ "no Plasmodium",
            is.na(sporozoites) ~ NA_character_ # no real need but just in case i messed up
        ),
        method = "dissection"
    ) |> # flag clearly for future potential debugging/differential treatment
    select(
        pot_unique, dpi, age, species, treatment, temp_mean, temp_range,
        replicate, pot, id, exp_age, label_pf, method
    ) # same column order as in surv_eip
# view(spz_diss)

surv_eip <- bind_rows(surv_eip, spz_diss) |>
    arrange(pot_unique, dpi) # combine and order properly. a row is a mosquito.

str(surv_eip)
surv_eip$species <- as.factor(surv_eip$species)
surv_eip$temp_mean <- as.factor(surv_eip$temp_mean)
surv_eip$temp_range <- as.factor(surv_eip$temp_range)
surv_eip$treatment <- as.factor(surv_eip$treatment)
surv_eip$pot <- as.factor(surv_eip$pot)
surv_eip$replicate <- as.factor(surv_eip$replicate)
surv_eip$label_pf <- as.factor(surv_eip$label_pf)
surv_eip$age <- as.integer(surv_eip$age)

# sum-to-zero contrasts
contrasts(surv_eip$species) <- contr.sum(nlevels(surv_eip$species))
contrasts(surv_eip$temp_mean) <- contr.sum(nlevels(surv_eip$temp_mean))
contrasts(surv_eip$temp_range) <- contr.sum(nlevels(surv_eip$temp_range))
contrasts(surv_eip$pot) <- contr.sum(nlevels(surv_eip$pot))


# Let's start with a GAMM
# not a GLM because we want flexibility in the shape of the logistic curve
# also, this allows us to test linarity of the predictors (edf values >1)

table(surv_eip$species, surv_eip$temp_mean, surv_eip$temp_range) # sample sizes are alright

m_gam <- gam(
    label_pf ~ s(age, by = pot) + # this fits a smooth curve for how age affects the probability of positivity, for each pot
        # label_pf ~ s(age, by = pot, k = 20) + # this fits a smooth curve for how age affects the probability of positivity, for each pot
        species * temp_mean * temp_range + # interaction terms fits a different intercept to each pot
        s(replicate, bs = "re"), # random effect for replicate
    family = binomial(link = "logit"), # logistic curve
    data = surv_eip,
    method = "REML" # see (Wood, 2011), also author of mgcv pkg
)
gam.check(m_gam) # diagnostic
# doubling k does not affect edf values
# from the docs (https://www.rdocumentation.org/packages/mgcv/versions/1.9-4/topics/gam.check and https://www.rdocumentation.org/packages/mgcv/versions/1.9-4/topics/choose.k)
# k-index: The further below 1 this is, the more likely it is that there is missed pattern left in the residuals (close to 1 is a good fit)
# but k-index/p-value analysis is not very robust for binary cases... Let's look at the k'/edf results instead.
# k' is the maximum curvature (df) allowed (default 10, but one is lost at fitting the smooth, so 9),
# edf (effective degrees of freedom, or 'smooth effects') is the curvature observed (edf=1 represents a linear relationship)
# edf < k' consistently, meaning the choice of k is alright
# random effect of replicate is very important
# plots reflect quirks of binary data. Residual signs match y-y_hat, y in {0, 1}, y_hat in(0, 1)
# (max residual for y=0 is 0, min residual for y=1 is 0)
# gap in (0, 1) reflects low frequency of high y_hat (see 4th plot).

summary(m_gam) # statistical inference
# most groups show nonlinear effect of age on eip cdf
# replicate is an important random effect

#              edf    p
# s(age):ci210   1    *
# s(age):ci216   1    ***
# s(age):ci270   4    ***
# s(age):ci276   2    .
# s(age):ki210   3    ***
# s(age):ki216   3    ***
# s(age):ki270   4    **
# s(age):ki276   3    **
# s(replicate)   7    ***

# plot results
gam_eip <- draw(m_gam, select = 1:8, ncol = 4)
gam_eip # this is cool

ggsave(plot = gam_eip, filename = "Figures/GAM_eip.png", width = 18, height = 12)
ggsave(plot = gam_eip, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/GAM_eip.png", width = 16, height = 12)

# conclusion:
# most cases show nonlinear growth
# no general evidence for decay (only some cases with sparse tail data, and 95%CI still includes an increase)
# general evidence for slowed down growth on tail (non linear relationship of age and p(infected) in the logit scale)
# temp treatments very similar for both species

# alright, let's fit a Bayesian model:
# EIP CDF should have a logistic shape
# GLM not valid (asymptote not 1, relationships not linear (GAMM result))
# replicate must be a random effect
# GAMM is very flexible and confirms shape makes sense, but limited interpretability.
# Given random effects, binary data (Bernouilli), and bounded parameters, lets do Bayesian


# we'll include the random effect only on 1 parameter of the final model...
# to decide which one, we'll try them all.
# let's define the brms formulae (bf) that we'll be trying out:

# 1. random effect on max prev
bform_upper <- bf( # builds a brms formula, to feed to brms later
    label_pf ~ inv_logit(logitupper) * inv_logit(exp(logslope) * (log(age) - logeip50)), # model formula. bounded to [0, 1] already. Product of max prev and standard logistic cdf
    logitupper ~ 0 + pot + (1 | replicate), # intercept for logit(asymptote)
    logslope ~ 0 + pot, # intercept for slope (around eip50, represents synchrony of sporogony. Modelled in logscale for numerical stability)
    logeip50 ~ 0 + pot, # intercept for eip50 (this is modelled in the log scale directly for numerical stability)
    nl = TRUE # nonlinear
)

# 2. random effect on eip50
bform_eip50 <- bf(
    label_pf ~ inv_logit(logitupper) * inv_logit(exp(logslope) * (log(age) - logeip50)),
    logitupper ~ 0 + pot,
    logslope ~ 0 + pot,
    logeip50 ~ 0 + pot + (1 | replicate),
    nl = TRUE
)

# 3. random effect on grotwth rate
bform_slope <- bf(
    label_pf ~ inv_logit(logitupper) * inv_logit(exp(logslope) * (log(age) - logeip50)),
    logitupper ~ 0 + pot,
    logslope ~ 0 + pot + (1 | replicate), # in log(age) scale, repre
    logeip50 ~ 0 + pot,
    nl = TRUE
)


# Also, we first need some priors
# we derive them from the GAMM results

# 1. create dataset for predictions
pot_lookup <- surv_eip |>
    distinct(pot, species, temp_mean, temp_range)

gamm_grid <- expand.grid(
    age = sort(unique(surv_eip$age)),
    pot = levels(surv_eip$pot)
) |>
    left_join(pot_lookup, by = "pot") |>
    mutate(replicate = factor(levels(surv_eip$replicate)[1])) # doesn't teally matter, won't be used. It just needs to exist as a colmun

# 2. add predictions according to the GAMM model
gamm_grid$fit <- predict(m_gam,
    newdata = gamm_grid, type = "response",
    exclude = "s(replicate)"
)

# 3. extract prior means from gamm result
group_ests <- gamm_grid |>
    group_by(pot) |>
    arrange(age, .by_group = TRUE) |> # eip50 and slope extractions below assume rows are in order of age
    summarise(
        asymptote = max(fit), # max prevalence
        eip50 = age[which(fit >= 0.5 * max(fit))[1]], # age at index 'first row where the predicted prev is > max/2'
        slope = {
            p_scaled <- pmin(fit / max(fit), 0.9999) # prevalence as a proportion of max
            logit_p <- qlogis(p_scaled) # ln(p/(1-p)) -> scaled prevalence onto logit scale, as modelled
            log_age <- log(age) # age onto ln scale, as modelled
            idx <- abs(age - eip50[1]) <= 1 # take those rows for which age is in eip50 ± 1 (in the real scale)
            coef(lm(logit_p[idx] ~ log_age[idx]))[2] # fit a linear regression (log-logit scale) and retrieve slope
        },
    )
# prev needs to be in logit scale, and age needs to be in log scale
# logit(p_scaled) normalises the asymptote, and straightens the curve.
# this is the logistic curve's growth rate parameter, not the curve itself.
# calculates sporogony rate relative to max prev and to relative passage of time
# it answers 'increase in odds of finding sporozoites per 1% age increase'


# visualise predictions one by one. priors come from these.
# They don't need to be perfect, just orientative (priors'll have wide variances)
toypot <- "ci216"
ggplot(gamm_grid |> filter(pot == toypot), aes(age, fit)) +
    geom_line() +
    geom_vline(xintercept = group_ests$eip50[group_ests$pot == toypot], linetype = "dashed") +
    labs(title = paste0(toypot, " GAM-fitted curve"))


# collapse to one shared, inflated prior per parameter
k <- 2.5 # inflate variances

# asymptote priors (normal in logit scale)
la <- qlogis(group_ests$asymptote) # predicted asymptotes to logit scale
mu_logitupper <- mean(la)
sd_logitupper <- k * sd(la)

# t priors (normal in log scale)
logeip50_pred <- log(group_ests$eip50) # predicted eip50 to log scale
mu_eip50 <- mean(logeip50_pred)
sd_eip50 <- k * sd(logeip50_pred)

# k priors (normal in log scale)
logslope_pred <- log(group_ests$slope) # predicted slope to log scale
mu_slope <- mean(logslope_pred)
sd_slope <- k * sd(logslope_pred)

# priors for fixed effects
# needs prior_string because values are stored in variables - brms handles prior specifications fragily with prior()
priors <- c(
    prior_string(sprintf("normal(%.2f,%.2f)", mu_logitupper, sd_logitupper), nlpar = "logitupper"),
    prior_string(sprintf("normal(%.2f,%.2f)", mu_eip50, sd_eip50), nlpar = "logeip50"),
    prior_string(sprintf("normal(%.2f,%.2f)", mu_slope, sd_slope), nlpar = "logslope")
)


# now create different sets of priors, for each case of random effect
# we'll compare them and choose the best one
priors_upper <- c(
    priors,
    prior_string("normal(0, 0.5)", class = "sd", nlpar = "logitupper", group = "replicate")
)

priors_eip <- c(
    priors,
    prior_string("normal(0, 0.5)", class = "sd", nlpar = "logeip50", group = "replicate")
)

priors_slope <- c(
    priors,
    prior_string("normal(0, 0.5)", class = "sd", nlpar = "logslope", group = "replicate")
)

# let's make sure these are mathematically sound and make sense
# the next few code blocks fit a cheaper version of the models
# that are fitted later, but without looking at the data,
# only sampling priors. This is to check:
# 1. that nothing's broken (models can run)
# 2. that the prior sampling space is okay (sampleable, not too narrow)

# g
m_prior_check_upper <- brm(bform_upper,
    data = surv_eip, family = bernoulli("identity"),
    prior = priors_upper, sample_prior = "only",
    chains = 2, iter = 2000
)
pp_check(m_prior_check_upper, type = "bars", ndraws = 50) # plot
draws <- as_draws_df(m_prior_check_upper)
cols <- grep("^b_logitupper_", names(draws), value = TRUE)
r_cols <- grep("^r_replicate__logitupper", names(draws), value = TRUE)
for (col in cols) {
    combined <- draws[[col]] + draws[[r_cols[1]]]
    cat(
        col, ": range on logit scale =", range(combined),
        "| range on probability scale =", range(plogis(combined)), "\n"
    )
} # samples between 0 and 1. good


# t
m_prior_check_eip <- brm(bform_eip50,
    data = surv_eip, family = bernoulli("identity"),
    prior = priors_eip, sample_prior = "only",
    chains = 2, iter = 2000
)
pp_check(m_prior_check_eip, type = "bars", ndraws = 50)
draws <- as_draws_df(m_prior_check_eip)
log_eip50_cols <- grep("^b_logeip50_", names(draws), value = TRUE)
r_cols <- grep("^r_replicate__logeip50", names(draws), value = TRUE)
for (col in log_eip50_cols) {
    combined <- draws[[col]] + draws[[r_cols[1]]]
    cat(
        col, ": range on log scale =", range(combined),
        "| range on eip50 scale =", range(exp(combined)), "\n"
    )
} # samples between <1 and stupidly high numbers. could be narrower, but good


# k
m_prior_check_slope <- brm(bform_slope,
    data = surv_eip, family = bernoulli("identity"),
    prior = priors_slope, sample_prior = "only",
    chains = 2, iter = 2000
)
pp_check(m_prior_check_slope, type = "bars", ndraws = 50)
draws <- as_draws_df(m_prior_check_slope)
cols <- grep("^b_logslope_", names(draws), value = TRUE)
r_cols <- grep("^r_replicate__logslope", names(draws), value = TRUE)
for (col in cols) {
    combined <- draws[[col]] + draws[[r_cols[1]]]
    cat(
        col, ": range on log scale =", range(combined),
        "| range on slope scale =", range(exp(combined)), "\n"
    )
} # all sampling is positive here too. also extremely wide but thats ok.

# all looks good.

# alright! Now the thing

# fit them - this steps takes ages, only do it once...
m_upper <- brm(bform_upper,
    data = surv_eip, family = bernoulli(link = "identity"),
    prior = priors_upper, chains = 4, cores = 4, iter = 4000,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
) # no warnings, clean convergence (r>4: 4 divergent transitions)

m_eip50 <- brm(bform_eip50,
    data = surv_eip, family = bernoulli(link = "identity"),
    prior = priors_eip, chains = 4, cores = 4, iter = 4000,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
) # clean convergence too (r>4: 11 divergent transitions)

m_slope <- brm(bform_slope,
    data = surv_eip, family = bernoulli(link = "identity"),
    prior = priors_slope, chains = 4, cores = 4, iter = 4000,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
) # same here (r>4: 3 divergent transitions, high rhat, low tail and bulk ess)

summary(m_upper)
summary(m_eip50)
summary(m_slope) # (r>4: >1 rhat, low tail and bulk ess)

# most important thing is that all Rhat are

loo(m_slope)
loo(m_eip50)
loo(m_upper)
# most important thing is that Pareto k estimates are between 0.5 and 0.7 (a handful of fails is alright, n=3010)

# not nested, and aic doesnt quite apply here cause random df <1. Compare with loo:
m_upper <- add_criterion(m_upper, "loo")
m_eip50 <- add_criterion(m_eip50, "loo")
m_slope <- add_criterion(m_slope, "loo")

comp <- loo_compare(m_upper, m_eip50, m_slope)
print(comp, digits = 2)

# loo comparison suggests that placing the random effect on upper makes most sense (indistinguishable from putting it on eip50)
# m_upper is the cleanest fit, and also makes most sense biologically.

bayesplot_theme_set(theme_default(base_size = 25)) # increase font size

pairs(m_slope, variable = c("b_logitupper_potci210", "b_logslope_potci210", "b_logeip50_potci210"))
pairs(m_eip50, variable = c("b_logitupper_potci210", "b_logslope_potci210", "b_logeip50_potci210"))
pairs(m_upper, variable = c("b_logitupper_potci210", "b_logslope_potci210", "b_logeip50_potci210"))

pairs(m_slope, variable = c("b_logitupper_potki210", "b_logslope_potki210", "b_logeip50_potki210"))
pairs(m_eip50, variable = c("b_logitupper_potki210", "b_logslope_potki210", "b_logeip50_potki210"))
pairs(m_upper, variable = c("b_logitupper_potki210", "b_logslope_potki210", "b_logeip50_potki210"))

pairs(m_slope, variable = c("b_logitupper_potki276", "b_logslope_potki276", "b_logeip50_potki276"))
pairs(m_eip50, variable = c("b_logitupper_potki276", "b_logslope_potki276", "b_logeip50_potki276"))
pairs(m_upper, variable = c("b_logitupper_potki276", "b_logslope_potki276", "b_logeip50_potki276"))

pairs(m_slope, variable = c("b_logitupper_potci276", "b_logslope_potci276", "b_logeip50_potci276"))
pairs(m_eip50, variable = c("b_logitupper_potci276", "b_logslope_potci276", "b_logeip50_potci276"))
pairs(m_upper, variable = c("b_logitupper_potci276", "b_logslope_potci276", "b_logeip50_potci276"))

# some of these posteriors seem highly entangled..

# save them all (helper for readability)
save_pair_plots <- function(model, model_name, pot) {
    vars <- c(
        paste0("b_logitupper_pot", pot),
        paste0("b_logslope_pot", pot),
        paste0("b_logeip50_pot", pot)
    )

    figname <- sprintf("pairs_%s_%s.png", model_name, pot)

    paths <- c(
        file.path("Figures", figname),
        file.path("/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics", figname)
    )

    plot <- pairs(model, variable = vars)

    for (p in paths) {
        png(p, width = 10, height = 10, units = "in", res = 300)
        print(plot)
        dev.off()
    }
}

models <- list(m_slope = m_slope, m_eip50 = m_eip50, m_upper = m_upper)
pots <- c("ci210", "ki210", "ki276", "ci276")

# takes ages
for (pot in pots) {
    for (mn in names(models)) {
        save_pair_plots(models[[mn]], mn, pot)
    }
}


# so, we specify our nonlinear logistic model
best_model <- m_upper

# check model fit
summary(best_model)
# RHAt should be 1 (defo not >1) - good
# BulkESS and Tail ESS in the thousands - good

# make predictions and explore effects as differences in predicted parameters
draws <- as_draws_df(best_model)
groups <- unique(surv_eip$pot)

get_param <- function(param, pot) { # convert params back to natural scale
    x <- draws[[paste0("b_", param, "_pot", pot)]]
    if (param == "logitupper") {
        plogis(x)
    } else if (param == "logeip50" | param == "logslope") {
        exp(x)
    } else {
        x
    }
}

# pull all 3 parameters for all 8 groups (for each draw, n=8000) into one long df, one row per draw x group
post <- map_dfr(groups, function(p) {
    slope <- get_param("logslope", p)
    eip50 <- get_param("logeip50", p)
    upper <- get_param("logitupper", p)
    tibble(
        pot = p,
        draw = seq_along(slope),
        asymptote = upper,
        eip50 = eip50,
        eip90_eip10 = eip50 * (9^(1 / slope) - 9^(-1 / slope)), # usual dispersion metric
        mass_80 = upper * 0.8, # so eip90-10 can be scaled to the amount of events happening in that time interval
        iqr = eip50 * (3^(1 / slope) - 3^(-1 / slope)), # iqr
        mass_50 = upper * 0.5 # so iqr can be scaled to the amount of events happening in that time interval
    )
})

# which comparisons are relevant
contrasts <- list(
    "species at 21±0" = c("ki210", "ci210"),
    "species at 21±6" = c("ki216", "ci216"),
    "species at 27±0" = c("ki270", "ci270"),
    "species at 27±6" = c("ki276", "ci276"),
    "Tmean for AC ±0" = c("ci210", "ci270"),
    "Tmean for AC ±6" = c("ci216", "ci276"),
    "Tmean for AG ±0" = c("ki210", "ki270"),
    "Tmean for AG ±6" = c("ki216", "ki276"),
    "Trange for AC 21" = c("ci210", "ci216"),
    "Trange for AC 27" = c("ci270", "ci276"),
    "Trange for AG 21" = c("ki210", "ki216"),
    "Trange for AG 27" = c("ki270", "ki276")
)

# helpre to build table - this summarises the parameter differences across the n=8000 draws
summarise_diff <- function(x) {
    n_total <- length(x)
    x <- x[is.finite(x)]
    n_dropped <- n_total - length(x) # drop cases where predicted diffs are Inf (if any)
    tibble(
        Median = median(x),
        CI_low = unname(quantile(x, .025)),
        CI_high = unname(quantile(x, .975)),
        pd_manual = max(mean(x > 0), mean(x < 0)), # probability of direction (proportion of draws on one side of 0)
        n_dropped = n_dropped # how many diffs were non-finite? (need to be discarded for calculating eip90-eip10)
    )
}

table <- map_dfr(names(contrasts), function(cname) { # for each relevant comparison (see contrasts)
    g1 <- contrasts[[cname]][1] # get the group names
    g2 <- contrasts[[cname]][2]
    d1 <- post %>% filter(pot == g1) # get their posteriors (see post)
    d2 <- post %>% filter(pot == g2)
    map_dfr(c("asymptote", "eip50", "eip90_eip10", "mass_80", "iqr", "mass_50"), function(p) { # for each of these colnames in post
        diff <- d1[[p]] - d2[[p]] # calculate their differences (in post, one row is one posterior draw)
        bind_cols(Comparison = cname, Parameter = p, summarise_diff(diff)) # summarise according to summarise_diff
    })
}) %>% # round and order rows
    mutate(across(c(Median, CI_low, CI_high, pd_manual), ~ round(., 3))) %>%
    arrange(Parameter, Comparison)

nrow(filter(table, n_dropped > 0)) # did any draw get discarded because the effect was Inf? (sign of instability in model)

# in case the above is >0
# table_finite <- table |>
#     filter(n_dropped == 0) |>
#     select(Comparison, Parameter, n_dropped)

# print(table_finite, n = Inf)

# check which pots contribute most to non-finite eip90_eip10 draws
post |>
    group_by(pot) |>
    summarise(
        n_nonfinite = sum(!is.finite(eip90_eip10)),
        pct = 100 * n_nonfinite / n()
    ) |>
    arrange(desc(pct))
# model with good fit prints all 0s here

# aaaaand plot things
################################################################################################################
set.seed(1984)

selected_colors <- c("#002fff", "#80b5ff", "#ff0000", "#ff87eb")

surv_exp <- filter(surv, exposed == "Exposed") # only exposed mosquitoes

# eip prediction grid
eip_grid <- expand.grid(
    age = sort(unique(surv_exp$age)), # we are gonna be doing log(age) later, but there's no 0s in age
    pot = levels(surv_eip$pot)
) |>
    left_join(pot_lookup, by = "pot") |> # add species, tmean, trange
    mutate(
        treatment = paste0(temp_mean, "±", temp_range, "°C"),
        row_id = row_number()
    )

# add EIP predictions from bayesian model
eip_preds <- posterior_epred(best_model, newdata = eip_grid, re_formula = NA)
B <- nrow(eip_preds)

# add to eip_grid (summarised)
eip_grid$predicted <- apply(eip_preds, 2, median)
eip_grid$lower <- apply(eip_preds, 2, quantile, probs = 0.25)
eip_grid$upper <- apply(eip_preds, 2, quantile, probs = 0.75)

# survival predictions from flexsurv model
#  first fit surv model (see 4_2_parametric_surv_exposed.r)
best_dist <- "genf" # from paraemtric survival modelling script
m1 <- flexsurvreg(Surv(age, dead) ~ species * temp_mean * temp_range, data = surv_exp, dist = best_dist) # from paraemtric survival modelling script

# dataset for predictions
surv_grid <- expand.grid(
    species = levels(surv_eip$species),
    temp_range = levels(surv_eip$temp_range),
    temp_mean = levels(surv_eip$temp_mean)
) |> left_join(pot_lookup, by = c("species", "temp_mean", "temp_range")) # add pot

# survival curve draws: one joint parameter draw per replicate -> one coherent, monotone curve
ages_grid <- sort(unique(eip_grid$age))
boot_pars <- normboot.flexsurvreg(m1, B = B, newdata = surv_grid) # bootstrap B parameter vectors for each pot value from the model (from their joint distribution so theyre coherent)

# make survival predictions (bootstrap)
S_t_array <- array(NA_real_, dim = c(B, length(ages_grid), nrow(surv_grid))) # 3 dims: bootstrap sample, age, pot
for (g in seq_len(nrow(surv_grid))) { # for each pot
    pars <- boot_pars[[g]] # get the bootstrapped genf parameters (n = B)
    for (b in seq_len(B)) { # for each parameter vector, populate the age dimension with genf draws
        S_t_array[b, , g] <- 1 - pgenf(ages_grid, # S(t) = 1-cdf, see chapter 2
            mu = pars[b, "mu"], sigma = pars[b, "sigma"],
            Q  = pars[b, "Q"],  P = pars[b, "P"]
        )
    }
}
# match eip and surv predicitons
group_idx <- match(eip_grid$pot, surv_grid$pot)
age_idx <- match(eip_grid$age, ages_grid)

# point predictions for plotting surv curves in the main figrue
S_t_median <- apply(S_t_array, c(2, 3), median) # median S(t) for each age*pot (dims 2, 3; now rows and cols)
eip_grid$S_t <- S_t_median[cbind(age_idx, group_idx)] # add to eip_grid

# for each row in eip_grid (pot*age, here treated as columns), record all survival preds (B rows)
S_t_draws <- sapply(seq_len(nrow(eip_grid)), function(i) S_t_array[, age_idx[i], group_idx[i]])

# now eip cdf predictions and survival estimates have the same structure
# 8k rows (bs draws) x 256 cols (pots x ages)
str(eip_preds)
str(S_t_draws)

# element-wise product of eip_preds * S_t_draws is the probability of being alive and infectious
area_draws <- eip_preds * S_t_draws
str(area_draws)
# for age*pot combination, and for every posterior draw from EIP and surv

# summarise area_draws (median ± iqr) and add to eip_grid
eip_grid$area_product <- apply(area_draws, 2, median) # median across 2nd dim (for each column)
eip_grid$area_product_lower <- apply(area_draws, 2, quantile, probs = 0.25)
eip_grid$area_product_upper <- apply(area_draws, 2, quantile, probs = 0.75)

area_df <- eip_grid # deserves renaming

# EIP posterior curves
################################################################################################################

eip_dists <- ggplot(eip_grid, aes(x = age, y = predicted, color = treatment, fill = treatment)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
    facet_wrap2(species ~ ., scales = "free") +
    scale_color_manual(values = selected_colors) +
    scale_fill_manual(values = selected_colors) +
    ylim(0, 0.8) +
    labs(
        x = "Mosquito age (days)",
        y = "Posterior prevalence ± IQR",
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
# slope itself means nothing, so instead it's used to calculate eip90-eip10

lookup <- surv_eip %>% # each model names the same thing a different way sfruayibuoehiabvhk
    distinct(species, temp_mean, temp_range, pot) %>%
    mutate(
        treatment = paste0(temp_mean, "±", temp_range)
    )

post_summary <- post %>% # posteriors from the eip model
    pivot_longer(c(asymptote, eip50, eip90_eip10), names_to = "Parameter", values_to = "value") %>%
    group_by(pot, Parameter) %>%
    summarise(
        Median = median(value),
        CI_low = quantile(value, .25),
        CI_high = quantile(value, .75),
        .groups = "drop"
    ) %>%
    left_join(lookup, by = "pot") %>%
    mutate(Parameter = factor(Parameter,
        levels = c("asymptote", "eip50", "eip90_eip10"),
        labels = c("Competence (%)", "EIP50 (d)", "EIP90-EIP10 (d)")
    ))

post_params <- ggplot(post_summary, aes(x = treatment, y = Median, color = treatment)) +
    geom_pointrange(aes(ymin = CI_low, ymax = CI_high), linewidth = 1, size = 2) +
    facet_grid2(species ~ Parameter, scales = "free", independent = "all") + # facet_wrap doesnt allow side titles :(
    labs(x = "Treatment", y = "Posterior estimate (median ± IQR)") +
    facetted_pos_scales(
        y = list(
            Parameter == "Competence (%)" ~ scale_y_continuous(limits = c(0, 1)),
            Parameter == "EIP50 (d)" ~ scale_y_continuous(limits = c(0, 50)),
            Parameter == "EIP90-EIP10 (d)" ~ scale_y_continuous(limits = c(0, 50))
        )
    ) +
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

# both curves and the area under their product - main figure 1
################################################################################################################

km1 <- get_km_data(surv_exp, surv_grid) # see functions.r

# unify levels so that the legends behave - ok maybe I shouldnt have defined the same thign a million times... no harm done.
treatment_levels <- sort(unique(as.character(sp_prev$treatment)))
sp_prev$treatment <- factor(sp_prev$treatment, levels = treatment_levels)
eip_grid$treatment <- factor(eip_grid$treatment, levels = treatment_levels)
km1$treatment <- factor(km1$treatment, levels = treatment_levels)
area_df$treatment <- factor(area_df$treatment, levels = treatment_levels)

surv_eips <- ggplot(combined_prev, aes(x = age, y = prevalence, colour = factor(treatment), shape = replicate)) +
    geom_point(size = 3) + # prevalence points (empirical)
    scale_shape_manual(values = 1:nlevels(sp_prev$replicate)) + # shaped by replicate
    geom_line( # parametric prevalence (eip CDF)
        data = eip_grid, aes(x = age, y = predicted, color = factor(treatment)),
        inherit.aes = FALSE, linewidth = 1
    ) +
    geom_line( # parametric survival (flexsurvreg)
        data = area_df,
        aes(x = age, y = S_t, colour = factor(treatment)),
        linewidth = 1, inherit.aes = FALSE
    ) +
    geom_step( # empirical KM: step function
        data = km1,
        aes(x = time, y = est, colour = factor(treatment)),
        linewidth = 1, inherit.aes = FALSE
    ) +
    geom_ribbon( # product of the two lines
        data = area_df,
        aes(x = age, ymin = 0, ymax = area_product, fill = factor(treatment)),
        inherit.aes = FALSE, alpha = 0.45
    ) +
    # geom_ribbon(data = area_df, # uncertainty around area - looks horrible
    #             aes(x = age, ymin = area_product_lower, ymax = area_product_upper, fill = factor(treatment)),
    #             inherit.aes = FALSE, alpha = 0.15, linetype = "dotted", color = NA) +
    facet_grid2(species ~ treatment,
        scales = "free",
        axes = "margins", remove_labels = "all",
    ) +
    scale_color_manual(values = selected_colors) +
    scale_fill_manual(values = selected_colors) +
    scale_y_continuous( # both axis have the same scale, but different names
        name = "Proportion alive at time t, S(t)",
        sec.axis = sec_axis(~., name = "Proportion infectious at time t, Prev(t)")
    ) +
    labs(
        x = "Age",
        colour = "Temperature",
        fill = "Temperature",
        shape = "Replicate"
    ) +
    guides( # format legend
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


# wee sensitivity check for the extrapolation of survival curves - how much of the area relies on this?
################################################################################################################

# list of the latest death observed per pot
last_death <- surv_exp %>%
    filter(dead == TRUE) %>%
    group_by(species, temp_mean, temp_range) %>%
    summarise(max_death_day = max(age))

# now compute how much of the area belongs to times > max_death_day for each pot
area_df_truncated <- area_df %>% # join area_df and last_death in a new df (to not pollute original)
    left_join(last_death, by = c("species", "temp_mean", "temp_range")) %>% # add last death col (max_death_day)
    group_by(species, treatment) %>% # for each pot
    summarise(
        prop_surv_extrapolated = sum(S_t[age > max_death_day]) / sum(S_t),
        area_total = sum(area_product), # total area
        area_extrapolated = sum(area_product[age > max_death_day]), # area mass of ages with extrapolated survival
        extr_proportion = area_extrapolated / area_total, # proportion of mass from extrapolated ages
        extr_diff = area_total - area_extrapolated, # absolute difference of total and extrapolated masses
    )

# print the useful part
print(area_df_truncated %>% select(species, treatment, area_total, extr_proportion))

# mean
weighted.mean(area_df_truncated$extr_proportion)

# infectious days - main figure 2
################################################################################################################

# quantify area under the product
groups_tbl <- eip_grid %>% # log row_ids of each pot (sp*treatment)
    group_by(species, treatment) %>%
    summarise(row_ids = list(row_id), .groups = "drop")

auc_summary <- pmap_dfr(groups_tbl, function(species, treatment, row_ids) { # using groups_tbl
    idx <- match(row_ids, eip_grid$row_id) # find relevant ids in eip_grid
    auc_draws <- rowSums(area_draws[, idx, drop = FALSE]) # for each row (draws), sum all columns (filtered by idx)
    tibble( # summarise for plotting
        species = species,
        treatment = treatment,
        AUC_median = median(auc_draws),
        AUC_low = quantile(auc_draws, 0.25),
        AUC_high = quantile(auc_draws, 0.75)
    )
})

print(auc_summary)

auc_plot <- ggplot(auc_summary, aes(x = treatment, y = AUC_median, color = treatment)) +
    geom_point(size = 8) +
    geom_errorbar(aes(ymin = AUC_low, ymax = AUC_high), width = 0.15, linewidth = 0.8) +
    facet_wrap(~species, scales = "free_x") +
    scale_color_manual(values = selected_colors) +
    labs(x = NULL, y = "Mean infectious days per mosquito ± IQR") +
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

# per-replicate plotting
################################################################################################################
# just plaiying around
# does per-replicate modelling improve the fit of the eip curves to the raw data?

obs_combos <- surv_eip %>% distinct(pot, replicate)

# prediction grid crossing age with the pot*replicate pairs
spaghetti_grid <- expand.grid(
    age = ages_grid,
    pot = levels(surv_eip$pot)
) |>
    left_join(pot_lookup, by = "pot") |>
    inner_join(obs_combos, by = "pot", relationship = "many-to-many") |> # one row per age*pot*replicate
    mutate(treatment = paste0(temp_mean, "±", temp_range, "°C"))

# rerun posterior draws with replicate random effect included (re_formulla = NULL)
preds_rep <- posterior_epred(best_model, newdata = spaghetti_grid, re_formula = NULL)
spaghetti_grid$predicted <- apply(preds_rep, 2, median)

# plot with thin lines per replicate, and a bold population-average curve
plot <- ggplot() +
    geom_point(
        data = sp_prev, aes(x = age, y = prevalence, color = treatment, size = sqrt(n), alpha = sqrt(n)),
        show.legend = TRUE
    ) +
    geom_line(
        data = spaghetti_grid,
        aes(x = age, y = predicted, group = interaction(pot, replicate), color = treatment),
        alpha = 0.35, linewidth = 0.4
    ) +
    facet_grid2(species ~ treatment,
        scales = "free",
        axes = "margins", remove_labels = "all",
    ) +
    geom_line(
        data = eip_grid,
        aes(x = age, y = predicted, color = treatment),
        linewidth = 1.3, linetype = "dashed"
    ) +
    scale_color_manual(values = selected_colors) +
    theme_minimal()

plot

# what's happening in (for example) ci270?
toy <- sp_prev |>
    filter(species == "An. coluzzii" & treatment == "27±0°C") |>
    select(age, replicate, prevalence, n) |>
    arrange(age, replicate)

view(toy)
