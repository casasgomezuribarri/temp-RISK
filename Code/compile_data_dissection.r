# process raw dissection csv file and generate operational datasets for analysis and plotting

##################################################################################################################################
# Environment
##################################################################################################################################
# load packages
packages <- c(
    "readxl",
    "readr",
    "dplyr",
    "tibble",
    "lubridate",
    "survival",
    "ggplot2",
    "knitr",
    "kableExtra",
    "webshot2",
    "janitor"
)
for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}

# a couple useful fucntions for setting up the right working directory:

# returns full path to parent folder of current script
whereami <- function() {
    this_file <- commandArgs() %>%
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file) == 0) {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

# set working directory to be the parent of the above
setwd_grandparent <- function() {
    this_script <- whereami()
    setwd(dirname(this_script))
    # print working directory
    print("Current working directory:")
    print(getwd())
}

# apply it
setwd_grandparent() # ../Survival (parent of the Code folder)

##################################################################################################################################
# dissection dataset
##################################################################################################################################

dissection <- read_csv("Data/dissection.csv")
# dissection <- dissection[, -1] # remove first column (row number)

# missing values are indicated by -1 in the following variables. change to NA
vars <- c("oocysts", "sporozoites_g1", "sporozoites_g2")
dissection <- dissection %>%
    mutate_at(vars(all_of(vars)), ~ ifelse(. == -1, NA, .))

# some format changes
dissection$cal_day <- dmy(dissection$cal_day)
dissection$midgut_dissection <- as.logical(dissection$midgut_dissection)
dissection$sporozoites_g1 <- as.logical(dissection$sporozoites_g1)
dissection$sporozoites_g2 <- as.logical(dissection$sporozoites_g2)


# new variable "sporozoites" is true if sporozoites present in either gland
dissection$sporozoites <- dissection$sporozoites_g1 | dissection$sporozoites_g2

# new variable pot_unique is a unique identifier for each pot (combination of replicate + pot (experimental group)
dissection$pot_unique <- paste(dissection$replicate, "_", dissection$pot, sep = "")

# visualise dataframe
# view(dissection) # there's a bug and NA's don't appear. Still useful for navigating the table
unique(dissection$sporozoites_g1) # but NAs are in the dataset
dissection[216, 8:12] # see? and they're perfectly carried over

# create new dataframe summarising dissection outcome per replicate:dpi:pot
##############################################################################

oocysts <- dissection %>%
    group_by(replicate, dpi, pot, pot_unique) %>% # makes a group for each unique combination of these variables (pot_unique is replicate+pot, but this ensures it is carried over)
    summarise( # for each combination of replicate, dpi, pot...
        n_midguts = sum(midgut_dissection, na.rm = TRUE), # number of midguts dissected
        prevalence = mean(oocysts != 0, na.rm = TRUE), # the mean of a boolean vector is the proportion of TRUE values
        oocyst_avg = mean(oocysts[oocysts != 0], na.rm = TRUE), # mean oocyst count
        oocyst_std = sd(oocysts[oocysts != 0], na.rm = TRUE),
        n_salivaryglands = sum(glands_dissection, na.rm = TRUE), # number of salivary glands dissected
        sporozoites = any(sporozoites, na.rm = TRUE), # any() performs a multiple OR operation - was there a mosquito with sporozoites?
        .groups = "drop"
    )

# add new variable :species
oocysts$species <- character(nrow(oocysts)) # initialise
# populate
for (i in 1:nrow(oocysts)) {
    if (substr(oocysts$pot[i], 1, 1) == "c") {
        oocysts$species[i] <- "An. coluzzii"
    } else if (substr(oocysts$pot[i], 1, 1) == "k") {
        oocysts$species[i] <- "An. gambiae"
    } else {
        print(paste("Error: species not recognised, see oocysts[", i, ", ]"))
    }
}

# add new variable :mean_temp
oocysts$mean_temp <- integer(nrow(oocysts)) # initialise
# populate
for (i in 1:nrow(oocysts)) {
    oocysts$mean_temp[i] <- substr(oocysts$pot[i], 3, 4) %>% as.numeric()
}

# add new variable :temp_range
oocysts$temp_range <- integer(nrow(oocysts)) # initialise
# populate
for (i in 1:nrow(oocysts)) {
    oocysts$temp_range[i] <- substr(oocysts$pot[i], 5, 5) %>% as.numeric()
}

# view(oocysts)

# create new dataframe summarising oocyst prevalence and counts per dpi:pot
##############################################################################

# basically oocysts but aggregating across replicates
prevalence <- dissection %>%
    group_by(dpi, pot) %>%
    summarise(
        n_midguts = sum(midgut_dissection, na.rm = TRUE), # number of midguts dissected
        mean_prevalence = mean(oocysts != 0, na.rm = TRUE), # the mean of a boolean vector is the proportion of TRUE values
        se_prevalence = sd(oocysts != 0, na.rm = TRUE) / sqrt(sum(!is.na(oocysts))), # se formula for a proportion, removing NAs
        mean_count = mean(oocysts[oocysts != 0], na.rm = TRUE), # mean oocyst count
        se_count = sd(oocysts[oocysts != 0], na.rm = TRUE) / sqrt(sum(oocysts != 0, na.rm = TRUE)), # se oocyst count, denominator is total non-zero, non-NA counts
        .groups = "drop"
    )

# add new variable :species
prevalence$species <- character(nrow(prevalence)) # initialise
# populate
for (i in 1:nrow(prevalence)) {
    if (substr(prevalence$pot[i], 1, 1) == "c") {
        prevalence$species[i] <- "An. coluzzii"
    } else if (substr(prevalence$pot[i], 1, 1) == "k") {
        prevalence$species[i] <- "An. gambiae"
    } else {
        print(paste("Error: species not recognised, see oocysts[", i, ", ]"))
    }
}

# add new variable :mean_temp
prevalence$mean_temp <- integer(nrow(prevalence)) # initialise
# populate
for (i in 1:nrow(prevalence)) {
    prevalence$mean_temp[i] <- substr(prevalence$pot[i], 3, 4) %>% as.numeric()
}

# add new variable :temp_range
prevalence$temp_range <- integer(nrow(prevalence)) # initialise
# populate
for (i in 1:nrow(prevalence)) {
    prevalence$temp_range[i] <- substr(prevalence$pot[i], 5, 5) %>% as.numeric()
}

# add new variable :treatment
prevalence$treatment <- character(nrow(prevalence)) # initialise
# populate
for (i in 1:nrow(prevalence)) {
    prevalence$treatment[i] <- paste0(prevalence$mean_temp[i], "±", prevalence$temp_range[i])
}
prevalence$treatment <- as.factor(prevalence$treatment)

# view(prevalence)

# create new dataframe summarising oocyst prevalence and counts per replicate:pot
##############################################################################

# to avoid zero-inflation, prevalence is calculated as
# the proportion of positive midguts dissected between t_o and t_s-1
# where t_o = day of first oocyst & t_s = day of first sporozoite

# before first oocyst, zeroes don't represent failed establishment
# some sporozoite positives had 0 oocysts, so by t_s some zeroes are false negs.

replicates <- dissection %>%
    group_by(replicate, pot) %>%
    mutate( # first add t_o and t_s for each replicte:pot
        t_o = dpi[which(oocysts != 0)][1],
        t_s = dpi[which(sporozoites == TRUE)][1]
    ) %>%
    filter(dpi >= t_o & dpi < coalesce(t_s, Inf)) %>% # drop rows where dpi < t_o | dpi >= t_s  (if t_s is NA, then no upper bound to dpi)
    summarise( # now calculate all metrics
        t_o = min(t_o),
        t_s = max(t_s),
        last_diss = max(dpi),
        n_midguts = sum(midgut_dissection, na.rm = TRUE),
        n_midguts_positive = sum(oocysts != 0, na.rm = TRUE),
        mean_prevalence = mean(oocysts != 0, na.rm = TRUE),
        se_prevalence = sd(oocysts != 0, na.rm = TRUE) / sqrt(sum(!is.na(oocysts))),
        mean_count = mean(oocysts[oocysts != 0], na.rm = TRUE),
        se_count = sd(oocysts[oocysts != 0], na.rm = TRUE) / sqrt(sum(oocysts != 0, na.rm = TRUE)),
        # .groups = "drop"
    )

# view(replicates)


# old code without dpi filtering:
# replicates2 <- dissection %>%
#     group_by(replicate, pot) %>%
#     summarise(
#         n_midguts = sum(midgut_dissection, na.rm = TRUE), # number of midguts dissected
#         n_midguts_positive = sum(oocysts != 0, na.rm = TRUE), # number of midguts with oocysts (oocysts == NA when midgut_dissection == FALSE)
#         mean_prevalence = mean(oocysts != 0, na.rm = TRUE), # the mean of a boolean vector is the proportion of TRUE values
#         se_prevalence = sd(oocysts != 0, na.rm = TRUE) / sqrt(sum(!is.na(oocysts))), # se formula for a proportion, removing NAs
#         mean_count = mean(oocysts[oocysts != 0], na.rm = TRUE), # mean oocyst count (only positive midguts)
#         se_count = sd(oocysts[oocysts != 0], na.rm = TRUE) / sqrt(sum(oocysts != 0, na.rm = TRUE)), # se oocyst count, denominator is total non-zero, non-NA counts
#         .groups = "drop"
#     )
# view(replicates2)

# add new variable :species
replicates$species <- character(nrow(replicates)) # initialise
# populate
for (i in 1:nrow(replicates)) {
    if (substr(replicates$pot[i], 1, 1) == "c") {
        replicates$species[i] <- "An. coluzzii"
    } else if (substr(replicates$pot[i], 1, 1) == "k") {
        replicates$species[i] <- "An. gambiae"
    } else {
        print(paste("Error: species not recognised, see oocysts[", i, ", ]"))
    }
}

# add new variable :mean_temp
replicates$mean_temp <- integer(nrow(replicates)) # initialise
# populate
for (i in 1:nrow(replicates)) {
    replicates$mean_temp[i] <- substr(replicates$pot[i], 3, 4) %>% as.numeric()
}

# add new variable :temp_range
replicates$temp_range <- integer(nrow(replicates)) # initialise
# populate
for (i in 1:nrow(replicates)) {
    replicates$temp_range[i] <- substr(replicates$pot[i], 5, 5) %>% as.numeric()
}

# add new variable :treatment
replicates$treatment <- character(nrow(replicates)) # initialise
# populate
for (i in 1:nrow(replicates)) {
    replicates$treatment[i] <- paste0(replicates$mean_temp[i], "±", replicates$temp_range[i])
}
replicates$treatment <- as.factor(replicates$treatment)

# view(replicates)

# create new dataframe summarising eip outcome per replicate:pot
##############################################################################

eip <- oocysts %>%
    group_by(replicate, pot, pot_unique) %>% # collapse over dpi
    summarise( # creates a new dataframe with one row per group. The columns will be as described below:
        species = species[1], # species is the same for all rows
        mean_temp = mean_temp[1], # mean_temp is the same for all rows
        temp_range = temp_range[1], # temp_range is the same for all rows
        treatment = paste(mean_temp[1], "\u00B1", temp_range[1]), # will be useful for plotting
        eip = dpi[which(sporozoites == TRUE)[1]], # the first dpi where sporozoites are present
        .groups = "drop"
    )

# add new variable max_span
eip$max_span <- numeric(nrow(eip))
for (i in 1:nrow(eip)) {
    pot <- eip$pot_unique[i] # get pot name
    eip$max_span[i] <- max(surv[surv$pot_unique == pot, "dpi"]) # maximum lifespan in survival data
    .groups <- "drop"
}

# generate survival estimates based on pot-wise KM
##############################################################################

# initialise new vars
eip$km_surv_at_eip <- numeric(nrow(eip)) # KM y-value at t_eip
eip$km_surv_at_eip_lower <- numeric(nrow(eip)) # KM y-value at t_eip + std.err
eip$km_surv_at_eip_upper <- numeric(nrow(eip)) # KM y-value at t_eip + std.err
eip$km_RMST_0toMax <- numeric(nrow(eip)) # integral of KM curve from t_0 to t_max
eip$km_RMST_0toMax_upper <- numeric(nrow(eip)) # integral of KM curve from t_0 to t_max + std.err
eip$km_RMST_0toMax_lower <- numeric(nrow(eip)) # integral of KM curve from t_0 to t_max - std.err
eip$km_RMST_EiptoMax <- numeric(nrow(eip)) # integral of KM curve from t_eip to t_max
eip$km_RMST_EiptoMax_upper <- numeric(nrow(eip)) # integral of KM curve from t_eip to t_max + std.err
eip$km_RMST_EiptoMax_lower <- numeric(nrow(eip)) # integral of KM curve from t_eip to t_max - std.err


for (i in 1:length(unique(eip$pot_unique))) { # for each pot represented in eip
    pot <- eip$pot_unique[i] # get pot name
    t_eip <- eip$eip[i] # get eip value
    print(pot)
    print(t_eip)

    # fit a km curve for the current pot
    kmfit_pot <- survfit(event ~ 1, data = surv[surv$pot_unique == pot, ])

    # extract survival estimate and lower & upper limits
    if (is.na(t_eip)) { # if sporozoites where not seen, all is NA
        eip$km_surv_at_eip[i] <- NA
        se <- NA
        eip$km_surv_at_eip_lower[i] <- NA
        eip$km_surv_at_eip_upper[i] <- NA
    } else { # if sporozoites were seen, estimate survival at that point
        eip$km_surv_at_eip[i] <- summary(kmfit_pot, times = t_eip)$surv
        se <- summary(kmfit_pot, times = t_eip)$std.err
        eip$km_surv_at_eip_lower[i] <- eip$km_surv_at_eip[i] + se
        eip$km_surv_at_eip_upper[i] <- eip$km_surv_at_eip[i] - se
    }


    # Extract full list of survival probabilities and times to calculate integrals
    surv_probs <- kmfit_pot$surv
    surv_times <- kmfit_pot$time
    surv_std_err <- kmfit_pot$std.err

    # Since the KM is a step function, the integral is just a sum of rectangular areas
    # Since step = 1, all rectangles have width = 1. The sum of areas = the sum of survival probabilities

    # Integral t_0 to t_max
    eip$km_RMST_0toMax[i] <- sum(surv_probs)
    eip$km_RMST_0toMax_upper[i] <- sum(surv_probs + surv_std_err)
    eip$km_RMST_0toMax_lower[i] <- sum(surv_probs - surv_std_err)

    # Integral t_eip to t_max (same as above but only when t >= t_eip)
    eip$km_RMST_EiptoMax[i] <- sum(surv_probs[surv_times >= t_eip])
    eip$km_RMST_EiptoMax_upper[i] <- sum(surv_probs[surv_times >= t_eip] + surv_std_err[surv_times >= t_eip])
    eip$km_RMST_EiptoMax_lower[i] <- sum(surv_probs[surv_times >= t_eip] - surv_std_err[surv_times >= t_eip])
}

# view(eip)

# format variables correctly
eip$replicate <- as.factor(eip$replicate)
eip$pot <- as.factor(eip$pot)
eip$pot_unique <- as.factor(eip$pot_unique)
eip$species <- as.factor(eip$species)
eip$mean_temp <- as.factor(eip$mean_temp)
eip$temp_range <- as.factor(eip$temp_range)


# create and save pretty table
pretty_table <- kable(eip, format = "html", escape = FALSE, table.attr = "class='table table-striped'") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
save_kable(pretty_table, file = "Figures/eip_data.html")
save_kable(pretty_table, file = "Figures/eip_data.jpeg")
