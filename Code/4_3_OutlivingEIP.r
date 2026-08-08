# apply TAUS and plot results
# what is the probability of surviving long enough to become infectious?

# author: Iván Casas

##################################################################################################################################
# Environment
##################################################################################################################################

source("Code/compile_data_survival.r") # sets the right env, process raw data, displays summary datasets, saves some useful descriptive stats

# install taus package
# install.packages("remotes")
# remotes::install_github("casasgomezuribarri/TAUS")
# install.packages(".", repos = NULL, type = "source")
library(TAUS) # might need to run the above if you're not me

# some packages
packages <- c(
  "RColorBrewer",
  "multcomp",
  "multcompView",
  "betareg",
  "dplyr",
  "tidyr",
  "effectsize",
  "survminer",
  "nortest",
  "ggsurvfit",
  "cramer",
  "ggtext",
  "ggnewscale"
)

for (i in packages) {
  if (!require(i, character.only = TRUE)) install.packages(i)
  library(i, character.only = TRUE)
}

#################################################################################
# load the data, initialise useful variables
#################################################################################

# only exposed are relevant
surv <- surv %>% filter(exposed == "Exposed")

# define a few useful things for TAUS
var <- "pot" # the main variable of interest
cats <- c("species", "temp_mean", "temp_range", "treatment") # other variables of interest
table(surv[[var]], useNA = "always") # check counts per level

time_var <- "age" # variable with time data
event_var <- "dead" # event variable
surv$event_factor <- factor(surv[[event_var]], levels = c(TRUE, FALSE), labels = c("Dead", "Censored")) # factor event variable for plotting
event_factor <- "event_factor"

max_time <- max(surv[[time_var]])
plot_names <- "phd_TAUS"

##################################################################################################################################
# TAUS
##################################################################################################################################

# compute conditional survival matrices
cond_surv <- cond_surv_mat(
  data = surv, # the dataset - a typical life table
  var = var, # variable of interest (string, must be categorical). Each level is a population for which the matrix is computed
  cats = cats, # other categorical variables that matter (vector of strings, all categorical variables).
  time_var = time_var, # the time variable (string)
  event_var = event_var, # the event variable (string)
  res = 1, # resolution of the conditional survival grids for t and tau
  conf_int_level = 0.95, # confidence intervals
  aggregate_by_cats = TRUE # if TRUE, will separate by each unique combination of values in var and cats (default FALSE)
)
# view(cond_surv$cond_surv)

# plot conditional matrices as heatmaps
cond_surv_heatmap <- cond_surv_plot(cond_surv, ncol_hm = 4)
# ggsave(plot = cond_surv_heatmap, paste0("Figures/", plot_names, "_heatmaps.png"), width = 6.5, height = 4.5, dpi = 400)

# these show the probability of surviving an age y=τ conditional on being age t=x
# green lines show last death observed in the experiment
##################################################################################################################################
#  plot transmission potential (Oτ at τ=EIP) for each group
##################################################################################################################################

selected_colors <- c("#0c36f6", "#5a8cd1", "#da0000", "#e976d6")

eip_summary <- eip |>
  group_by(pot) |>
  summarise(
    mean_eip = mean(eip, na.rm = TRUE),
    se_eip = sd(eip, na.rm = TRUE) / sqrt(sum(!is.na(eip))), # se = sd / sqrt(n)
    .groups = "drop"
  ) |>
  mutate(treatment = paste0(substr(pot, 3, 4), "±", substr(pot, 5, 5)))


# get var_values
var_values <- unique(cond_surv$cond_surv$unique_label) # list
# get tau_values
tau_values <- floor(eip_summary$mean_eip[match(substr(var_values, 1, 5), eip_summary$pot)])

target_groups <- tibble(unique_label = var_values, tau = tau_values)


# filter the dataset
eip_tau <- cond_surv$cond_surv %>%
  semi_join(target_groups, by = c("unique_label", "tau"))
eip_tau$treatment_label <- eip_tau$treatment

# stats
stat_results <- pairwise_test(cond_surv,
  var_values = var_values,
  tau_values = tau_values
)
stat_results
# how to interpret the results:

#    O_tau_n = probability of outliving tau_n for a randomly selected individual from group_n
#    effect_size = O_tau_1 - O_tau_2
#      an individual from group_1 has effect_size more chance of outliving tau_1 than an individual from group_2 of outliving tau_2
#    effect_ratio = O_tau_1 / O_tau_2
#      an individual from group_1 has effect_ratio times the chance of outliving tau_1 than an individual from group_2 of outliving tau_2


#    ks statistic:
#    - 0: both distributions are identical
#    - 1: both distributions are completely different
#    ks statistic is analytically calculated as follows:
#    1. Betas are fitted to both O_tau (±95%CI) values, estimating alpha and beta parameters
#    2. ks = max(F1(x) - F2(x)) for all x in [0, 1], where F1 and F2 are the CDFs of the two Betas


#    p-value is likelihood of H0 being true (both distributions are similar)
#    calculated as follows
#    1. A null Beta is estimated by averaging the alpha and beta of the two Betas fitted in the calculation of the ks statistic
#    2. Pairs of values are sampled from the null Beta
#    3. the ks distance from those two sets of values is estimated and stored
#    4. Steps 3-4 are repeated B (=5000) times
#    5. p-value is the proportion of ks distances that are greater than the analytic ks.

#    p-value answers:
#    how often is the KS statistic from the null pairs greater than or equal to the deterministic one?


# default behaviour is to compare all possible pairs (Tukey-like)
# see only the relevant ones

# effect of species
stat_results_sp <- stat_results |>
  filter(substr(group_1, 2, 5) == substr(group_2, 2, 5)) |> # only include if everything is the same except species
  mutate(
    pot_1 = substr(group_1, 1, 5),
    pot_2 = substr(group_2, 1, 5)
  )
stat_results_sp

# effect of tm
stat_results_tm <- stat_results |>
  filter(paste0(substr(group_1, 1, 1), substr(group_1, 5, 5)) == paste0(substr(group_2, 1, 1), substr(group_2, 5, 5))) |> # only include if everything is the same except tm
  mutate(
    pot_1 = substr(group_1, 1, 5),
    pot_2 = substr(group_2, 1, 5)
  )
stat_results_tm

# effect of tr
stat_results_tr <- stat_results |>
  filter(paste0(substr(group_1, 1, 1), substr(group_1, 4, 4)) == paste0(substr(group_2, 1, 1), substr(group_2, 4, 4))) |> # only include if everything is the same except tr
  mutate(
    pot_1 = substr(group_1, 1, 5),
    pot_2 = substr(group_2, 1, 5)
  )
stat_results_tr


# manual annotation of significant differences
sig_labels <- data.frame(
  treatment_label = unique(eip_tau$treatment_label),
  x_pos = 1.5, # x position in plot
  y_pos = 0.95 # y position in plot
) |>
  mutate(treatment = paste0(
    substr(treatment_label, 1, 2),
    substr(treatment_label, 4, 4)
  )) |>
  rowwise() |>
  mutate(
    p_value = stat_results_sp[substr(stat_results_sp$pot_1, 3, 5) == treatment, ]$p_value,
    label = ifelse(p_value < 0.001, "***",
      ifelse(p_value < 0.01, "**",
        ifelse(p_value < 0.05, "*",
          ""
        )
      )
    )
  ) |>
  ungroup()
sig_labels

otau_sp <- ggplot(eip_tau, aes(x = species, y = O_tau, color = treatment, shape = species)) +
  geom_point(position = position_dodge(width = 0.5), size = 8) +
  geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.4, position = position_dodge(width = 0.5), linewidth = 1) +
  labs(x = "", y = "Probability of outliving EIP", title = "") +
  ylim(0, 1) +
  scale_color_manual(values = selected_colors) +
  theme_minimal(base_size = 18) +
  geom_text(
    data = sig_labels, aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE, size = 12
  ) +
  facet_wrap(~treatment_label, scales = "free", ncol = 4) +
  guides(
    color = guide_legend(ncol = 2, title = "Treatment", title.position = "top"),
    shape = guide_legend(ncol = 1, title = "Species", title.position = "top")
  ) +
  theme(
    plot.title = element_text(size = 39, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(size = 35),
    axis.text = element_text(size = 26, angle = 25, hjust = 1),
    strip.text = element_text(size = 28, face = "bold"),
    # strip.background = element_rect(fill = "#e8e7ff", colour = "#8c8cff", linewidth = 2),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 36),
    legend.position = "bottom",
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    panel.grid.major = element_line(color = "#eaeaea"),
    panel.grid.minor = element_blank()
  )

otau_sp

# save plot - this is the supplementary figure
ggsave(plot = otau_sp, paste0("Figures/", plot_names, "_P(T>tau)otau_sp.png"), width = 23, height = 12)
ggsave(plot = otau_sp, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/otau_sp.png", width = 16, height = 12)

str(eip_tau)

# no show all temps by species
otaus <- ggplot(eip_tau, aes(x = treatment, y = O_tau, color = treatment, shape = species)) +
  geom_point(position = position_dodge(width = 0.5), size = 8) +
  geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.4, position = position_dodge(width = 0.5), linewidth = 1) +
  labs(x = "", y = "Probability of outliving EIP", title = "") +
  ylim(0, 1) +
  scale_color_manual(values = selected_colors) +
  scale_x_discrete(labels = function(x) substr(x, 1, 4)) + # quick fix to remove the degree celsius symbol
  theme_minimal() +
  facet_wrap(~species, scales = "free") +
  guides(
    color = guide_legend(ncol = 2, title = "Treatment", title.position = "top"),
    shape = guide_legend(ncol = 1, title = "Species", title.position = "top")
  ) +
  theme(
    # plot.title = element_text(size = 39, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(size = 35),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 35, face = "bold"),
    # strip.background = element_rect(fill = "#e8e7ff", colour = "#8c8cff", linewidth = 2),
    legend.position = "none",
    # legend.text = element_text(size = 28),
    # legend.title = element_text(size = 36),
    # legend.position = "right",
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    panel.grid.major = element_line(color = "#eaeaea"),
    panel.grid.minor = element_blank()
  )
otaus

# save plot
ggsave(plot = otaus, paste0("Figures/", plot_names, "_P(T>tau)otau.png"), width = 16, height = 12)
ggsave(plot = otaus, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/otau.png", width = 10, height = 7.5)

# extract legend for plotting
otaus_legend <- otaus + theme(
  legend.position = "right",
  legend.text = element_text(size = 15), # dont wanna change the theme in the plot
  legend.title = element_text(size = 20)
) # but this legend text is tiny!
legend_only <- cowplot::get_legend(otaus_legend)
cowplot::save_plot(
  "/Users/ivancasas/GitHub/Thesis/Chapters/04_RISK/pics/colormap.png",
  legend_only,
  base_width = 4,
  base_height = 3
)
