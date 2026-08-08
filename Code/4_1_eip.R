# plot eip25
# glm (nb) on eip25 ~ species + mean temp * temp range

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
    "emmeans",
    "cowplot",
    "nlme",
    "flexsurv",
    "dplyr",
    "MASS",
    "tidyverse",
    "betareg",
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
# EIP by treatment (first sporozoite observed)
##################################################################################################################################
selected_colors <- c("#002fff", "#80b5ff", "#ff0000", "#ff87eb")

# view(eip)
eip$treatment <- as.factor(eip$treatment)

eip_summary <- eip |>
    group_by(pot) |>
    summarise(
        mean_eip = mean(eip, na.rm = TRUE),
        se_eip = sd(eip, na.rm = TRUE) / sqrt(sum(!is.na(eip))), # se = sd / sqrt(n)
        .groups = "drop"
    ) |> # for pretty plot
    mutate(
        treatment = paste0(substr(pot, 3, 4), "±", substr(pot, 5, 5)),
        species = case_when(
            str_sub(pot, 1, 1) == "k" ~ "An. gambiae",
            str_sub(pot, 1, 1) == "c" ~ "An. coluzzii"
        ),
        plot_label = paste0(species, ", ", treatment)
    )

eips <- ggplot(eip_summary, aes(x = treatment, y = mean_eip, color = treatment, shape = species)) +
    scale_y_continuous(limits = c(0, 30)) +
    geom_point(size = 10) +
    geom_errorbar(
        aes(ymin = mean_eip - se_eip, ymax = mean_eip + se_eip),
        width = 0.6
    ) +
    # geom_vline(xintercept = 4.5, color = "black", alpha = 0.3) + # add vertical line to separate species
    # geom_vline(xintercept = 0.5, color = "black", alpha = 0.3) + # add vertical line to separate species
    scale_color_manual(values = selected_colors) +
    scale_x_discrete() +
    theme_minimal() +
    facet_wrap(~species, scales = "free", ncol = 2) +
    labs(
        # title = "EIP across groups", # no title in publications
        x = "Temperature",
        y = "EIP (days until first sporozoite)",
        color = "Treatment",
        shape = "Species"
    ) +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 35),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(colour = "black", linewidth = 0.6),
        plot.margin = margin(10, 10, 10, 10),
        strip.text = element_text(size = 35, face = "bold"),
        legend.position = "none",
        # plot.title = element_text(size = 35, hjust = 0.5),
        # legend.title = element_text(size = 30),
        # legend.text = element_text(size = 24),
        # legend.key.size = unit(1.5, "cm")
    )

eips
ggsave(filename = "Figures/eips.png", plot = eips, width = 16, height = 12)
ggsave(filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/eips.png", plot = eips, width = 12, height = 9)

# stats (GLM)
##################################################################################################################################

str(eip)

contrasts(eip$species) <- contr.sum(nlevels(eip$species))
contrasts(eip$mean_temp) <- contr.sum(nlevels(eip$mean_temp))
contrasts(eip$temp_range) <- contr.sum(nlevels(eip$temp_range))

eip_model <- glm.nb(eip ~ species + mean_temp * temp_range, data = eip)


Anova(eip_model, type = "III") # for pvalues
summary(eip_model) # for coefs
