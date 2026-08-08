# plot prevalences and intensities against time (oocysts abd sporozoites)
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
    "merTools",
    "emmeans",
    "cowplot",
    "nlme",
    "flexsurv",
    "dplyr",
    "tidyverse",
    "betareg",
    "GGally",
    "ggh4x",
    "ggh4x",
    "lme4",
    "car"
)
for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}

# EIP parametric predictions
##################################################################################################################################
sp_prev <- surv |> # from 4_4
    filter(exposed == "Exposed") |>
    group_by(pot_unique, age, species, treatment, temp_mean, temp_range, replicate) |>
    summarise(
        mean_prevalence = mean(label_pf == "with Plasmodium", na.rm = TRUE), # proportion of positives among those qpcr'ed
        n_qpcr = sum(!is.na(label_pf)),
        mean_count = mean(pf_copies[pf_copies != 0], na.rm = TRUE) # mean dna count among all positives
    ) |>
    ungroup()

# grid
treatment_lookup <- surv %>%
    distinct(treatment, temp_mean, temp_range)
newdata <- expand.grid(
    age = seq(0, 30, length.out = 100),
    species = levels(surv$species),
    treatment = levels(surv$treatment)
) %>%
    left_join(treatment_lookup, by = "treatment")

# eip predictions
newdata$predicted <- predict(best_model,
    newdata = newdata,
    type = "response", re.form = NA
)

# parametric survival predictions
##################################################################################################################################

surv_exp <- filter(surv, exposed == "Exposed") # only with exposed mosquitoes
best_dist <- "genf" # from 4_1
m1 <- flexsurvreg(event ~ species * temp_mean * temp_range, data = surv_exp, dist = best_dist) # from 4_1

# grid
nd1 <- bind_rows(
    expand.grid(
        species = c("An. gambiae", "An. coluzzii"),
        temp_range = as.factor(c(0, 6)),
        temp_mean = as.factor(c(21, 27)),
        pfstatus_thr1 = c("Control", "Exposed")
    )
)

# predictions
pred1 <- summary(m1, newdata = nd1, type = "survival", ci = TRUE, tidy = TRUE)
pred1 <- pred1 %>%
    mutate(treatment = paste0(temp_mean, "±", temp_range, "°C")) |>
    rename(age = time) |> # align for later
    group_by(species, temp_mean, temp_range) %>%
    mutate(S_t = est) %>%
    ungroup()

# Infectious days
##################################################################################################################################
area_df <- newdata %>%
    inner_join(pred1, by = c("age", "species", "treatment", "temp_mean", "temp_range")) %>%
    mutate(area_product = predicted * S_t)

nrow(newdata)
nrow(pred1)
nrow(area_df)



# sensitivity analysis for of S(t) and Prev(t) on the result of the area value
