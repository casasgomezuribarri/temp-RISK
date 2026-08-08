# This script takes the raw excel files and compiles datasets that are ready for survival analysis

# author: Ivan Casas

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
  "scales",
  "htmltools",
  "kableExtra",
  "webshot2",
  "janitor",
  "tidyverse"
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
# survival dataset (surv)
##################################################################################################################################

# The data is in the first sheet of an excel file, so let's make it a standalone .csv first
rawdatapath_surv <- "Data/survival.xlsx"
data_surv <- read_excel(rawdatapath_surv, sheet = 1)
write.csv(data_surv, "Data/survival.csv")

# load data
datapath_surv <- "Data/survival.csv"
raw_surv <- read_csv(datapath_surv, show_col_types = FALSE)

# Some format changes
raw_surv$date <- ymd(raw_surv$cal_day)
raw_surv$dead <- as.integer(raw_surv$dead)
raw_surv$replicate <- as.factor(raw_surv$replicate)
raw_surv$pot <- as.factor(raw_surv$pot)
raw_surv$menu <- as.factor(raw_surv$menu)
raw_surv$alive <- as.integer(raw_surv$alive)
# str(raw_surv)

# initialize an empty list to store the life-table dataframes
surv_list <- list()

# convert to life table
for (i in 1:nrow(raw_surv)) {
  row <- raw_surv[i, ]

  # exposure
  exp <- switch(substr(row$pot, 2, 2),
    "i" = "Exposed",
    "c" = "Control",
    "o" = "Control"
  )

  # temp_range
  t_r <- as.integer(substr(row$pot, 5, 5))

  # species
  sp <- switch(substr(row$pot, 1, 1),
    "k" = "An. gambiae",
    "c" = "An. coluzzii",
  )

  # temp_mean
  t_m <- as.integer(substr(row$pot, 3, 4))

  # pot unique ID
  pot_unique <- as.factor(paste(row$replicate, "_", row$pot, sep = "")) # Pot value

  # dataframe for dead mosquitoes with right censored = FALSE
  dead_df <- data.frame(
    dpi = rep(row$dpi, row$dead),
    date = rep(row$date, row$dead),
    dead = rep(TRUE, row$dead),
    exposed = rep(exp, row$dead),
    species = rep(sp, row$dead),
    temp_mean = rep(t_m, row$dead),
    temp_range = rep(t_r, row$dead),
    treatment = rep(paste0(t_m, "±", t_r, "°C"), row$dead), # will be useful for plotting
    pot = rep(row$pot, row$dead),
    replicate = rep(row$replicate, row$dead),
    pot_unique = rep(pot_unique, row$dead)
  )

  # dataframe for alive mosquitoes with right censored = TRUE
  alive_df <- data.frame(
    dpi = rep(row$dpi, row$alive),
    date = rep(row$date, row$alive),
    dead = rep(FALSE, row$alive),
    exposed = rep(exp, row$alive),
    species = rep(sp, row$alive),
    temp_mean = rep(t_m, row$alive),
    temp_range = rep(t_r, row$alive),
    treatment = rep(paste0(t_m, "±", t_r, "°C"), row$alive), # will be useful for plotting
    pot = rep(row$pot, row$alive),
    replicate = rep(row$replicate, row$alive),
    pot_unique = rep(pot_unique, row$alive)
  )

  # Combine dead and alive dataframes into a single dataframe
  combined_df <- rbind(dead_df, alive_df)

  # Add the combined dataframe to the list
  surv_list[[i]] <- combined_df
}

surv <- do.call(rbind, surv_list) # format as a single dtaframe

# add replicate and pot_unique ids
surv$rep_int <- as.integer(factor(surv$replicate)) # to use replicate as frailty term
surv$pot_int <- as.integer(factor(surv$pot_unique)) # to use pot as frailty term

# add sample id (for maching with qPCR data)
surv <- surv %>%
  group_by(pot_unique, dead, dpi) %>%
  mutate(id = row_number()) %>% # first a mosquito ID within each pot:dpi:dead group
  ungroup() %>%
  mutate(sample = paste0( # now sample id code
    replicate,
    "_",
    toupper(pot), # pot in uppercase
    "_",
    sprintf("%02d", dpi), # age with 2 digits
    ".",
    id,
    ".",
    case_when(
      dead == TRUE ~ "x",
      dead == FALSE ~ "0"
    )
  )) %>% # let's calculate age = dpi + mean age at exposure
  group_by(replicate, species) %>%
  mutate(
    exp_age = case_when(
      species == "An. gambiae" & replicate == "r00" ~ 4, # AG r0: 4
      species == "An. gambiae" & replicate == "r01" ~ 4, # AG r1: 4
      species == "An. gambiae" & replicate == "r02" ~ 2, # AG r2: 2
      species == "An. gambiae" & replicate == "r03" ~ 2, # AG r3: 2
      species == "An. gambiae" & replicate == "r04" ~ 4, # AG r4: 4
      species == "An. gambiae" & replicate == "r05" ~ 3, # AG r5: 3
      species == "An. gambiae" & replicate == "r06" ~ 4, # AG r6: 4
      species == "An. gambiae" & replicate == "r07" ~ 5, # AG r7: 5
      species == "An. gambiae" & replicate == "r08" ~ 4, # AG r8: 4
      species == "An. gambiae" & replicate == "r09" ~ 4, # AG r9: 4
      species == "An. coluzzii" & replicate == "r00" ~ 3, # AC r0: 3
      species == "An. coluzzii" & replicate == "r01" ~ 0, # AC r1: 0
      species == "An. coluzzii" & replicate == "r02" ~ 2, # AC r2: 2
      species == "An. coluzzii" & replicate == "r03" ~ 2, # AC r3: 2
      species == "An. coluzzii" & replicate == "r04" ~ 4, # AC r4: 4
      species == "An. coluzzii" & replicate == "r05" ~ 1, # AC r5: 1
      species == "An. coluzzii" & replicate == "r06" ~ 3, # AC r6: 3
      species == "An. coluzzii" & replicate == "r07" ~ 5, # AC r7: 5
      species == "An. coluzzii" & replicate == "r08" ~ 4, # AC r8: 4
      species == "An. coluzzii" & replicate == "r09" ~ 4, # AC r9: 4
    )
  ) %>%
  ungroup() %>%
  mutate(
    age = dpi + exp_age
  )

surv$event <- with(surv, Surv(age, dead)) # event times (formality for dealing with right censored data)


##################################################################################################################################
# add oocyst prevalence and count
##################################################################################################################################
# this was measured only for some pots, so it is estimated based on treatment

# process dissection data
source("Code/compile_data_dissection.r") # process dissection data, shows all different operational datasets created

# dissectoin dataset has the raw dissectoin data. aggregate by pot_unique
oocysts_treatment <- dissection %>%
  group_by(pot_unique) %>%
  summarise(
    prevalence_oo = mean(oocysts != 0, na.rm = TRUE), # the mean of a boolean vector is the proportion of TRUE values
    prevalence_oo_se = sd(oocysts != 0, na.rm = TRUE) / sqrt(sum(!is.na(oocysts))), # se formula for a proportion, removing NAs
    count_oo = mean(oocysts[oocysts != 0], na.rm = TRUE), # mean oocyst count
    count_oo_se = sd(oocysts[oocysts != 0], na.rm = TRUE) / sqrt(sum(oocysts != 0 & !is.na(oocysts))), # se formula for a mean, removing NAs
    midgut_n = sum(midgut_dissection)
  ) %>% # now add the pot column again for the predicted ones
  left_join(distinct(dissection, pot_unique, pot), by = "pot_unique") %>%
  group_by(pot) %>%
  mutate(across(
    c(prevalence_oo, prevalence_oo_se, count_oo, count_oo_se),
    ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
  )) %>%
  ungroup()

# view(oocysts_treatment)

# add oocyst prevalences and counts to surv
surv <- surv %>%
  group_by(pot_unique) %>%
  left_join(select(oocysts_treatment, -pot), by = "pot_unique") %>%
  ungroup()

# started dissecting at r4. Pots that weren't dissected get their prevalences estimated
surv <- surv %>%
  group_by(pot) %>%
  mutate(
    prevalence_oo = ifelse(is.na(prevalence_oo), mean(prevalence_oo, na.rm = TRUE), prevalence_oo),
    prevalence_oo_se = ifelse(is.na(prevalence_oo_se), mean(prevalence_oo_se, na.rm = TRUE), prevalence_oo_se),
    count_oo = ifelse(is.na(count_oo), mean(count_oo, na.rm = TRUE), count_oo),
    count_oo_se = ifelse(is.na(count_oo_se), mean(count_oo_se, na.rm = TRUE), count_oo_se),
    midgut_n = ifelse(is.na(midgut_n), 0, midgut_n) # na here is no dissections so 0
  ) %>%
  ungroup()

##################################################################################################################################
# add sporozoite prevalence and intensity
##################################################################################################################################

# import qpcr data (sporozoites)
qpcr_data <- read_csv("Data/qpcrs_clean.csv", show_col_types = FALSE)
# str(qpcr_data)

qpcr_data %>%
  count(unique_ID) %>%
  filter(n > 1) # no duplicates on the qpcr id
surv %>%
  count(sample) %>%
  filter(n > 1) # no duplicates on the surv id


# add pf_label and pf_copies to surv
surv <- surv %>%
  left_join(qpcr_data %>% select(unique_ID, label_pf, pf_copies),
    # by = c("replicate", "sample")
    by = join_by("sample" == "unique_ID"),
  ) %>%
  mutate(label_pf = as.factor(label_pf))
table(surv$label_pf) # there should be 236 with Plasmodium and 1918 without

# hm some clean ones aren't matched
qpcr_data$unique_ID[!(qpcr_data$unique_ID %in% surv$sample)] # ah, they're from that odd CO276 pot, that's ok, none of those was actually exposed

qpcr_data$unique_ID[
  !(qpcr_data$unique_ID %in% surv$sample) & !grepl("CO", qpcr_data$unique_ID, fixed = TRUE)
] # there's 9 left
# purrrrrfact

# calculate prevalence and intensities of sporozoites by pot
surv <- surv %>%
  group_by(pot_unique) %>%
  mutate(
    prevalence_sp = mean(label_pf == "with Plasmodium", na.rm = TRUE), # proportion of positives among those qpcr'ed
    n_qpcr = sum(!is.na(label_pf)), # this will be useful for weighted means later on
    load_sp = mean(pf_copies[label_pf == "with Plasmodium"], na.rm = TRUE), # mean copies in those that tested positive
    n_pos_qpcr = sum(label_pf == "with Plasmodium", na.rm = TRUE) # number of positive qpcr results
  ) %>%
  ungroup()

##################################################################################################################################
# print and save dataset
##################################################################################################################################

# view(surv)
write.csv(surv, file = "Data/surv.csv", row.names = FALSE)

##################################################################################################################################
# descriptive stats
##################################################################################################################################

cm_pots <- ggplot(surv, aes(x = pot)) +
  geom_bar() +
  labs(title = "Count Map of Pots", x = "Pot", y = "Count") +
  theme_minimal() +
  theme(
    # plot.title = element_text(size = 30)
    # , axis.text.x = element_text(size = 24)
  )

ggsave("Figures/count_map_pots.png", plot = cm_pots, width = 8, height = 6)

pf_pots <- ggplot(surv, aes(x = pot, y = prevalence_oo)) +
  geom_boxplot() +
  labs(title = "Prevalence by Pot", x = "Pot", y = "Prevalence") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24)
  )
ggsave("Figures/prevalence_by_pots.png", plot = pf_pots, width = 8, height = 6)


cm_tempmean <- ggplot(surv, aes(x = temp_mean)) +
  geom_bar() +
  labs(title = "Count Map of Mean temperature", x = "Mean tempreature", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24)
  )

ggsave("Figures/count_map_tempmean.png", plot = cm_tempmean, width = 8, height = 6)

cm_temprange <- ggplot(surv, aes(x = temp_range)) +
  geom_bar() +
  labs(title = "Count Map of Temperature range", x = "Tempreature range", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24)
  )

ggsave("Figures/count_map_temprange.png", plot = cm_temprange, width = 8, height = 6)


cm_species <- ggplot(surv, aes(x = species)) +
  geom_bar() +
  labs(title = "Count Map of Species", x = "Species", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24)
  )

ggsave("Figures/count_map_species.png", plot = cm_species, width = 8, height = 6)


cm_exposure <- ggplot(surv, aes(x = exposed)) +
  geom_bar() +
  labs(title = "Count Map of Exposure", x = "Exposure", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24)
  )

ggsave("Figures/count_map_exposure.png", plot = cm_exposure, width = 8, height = 6)

################################################################################################
# pretty table to show sample sizes by pot and replicate
################################################################################################
surv_summary <- surv %>%
  group_by(replicate) %>%
  summarise(
    infection_date = first(date),
    length_days = as.character(max(dpi)),
    kc210 = sum(pot == "kc210"),
    ki210 = sum(pot == "ki210"),
    cc210 = sum(pot == "cc210"),
    ci210 = sum(pot == "ci210"),
    kc216 = sum(pot == "kc216"),
    ki216 = sum(pot == "ki216"),
    cc216 = sum(pot == "cc216"),
    ci216 = sum(pot == "ci216"),
    kc270 = sum(pot == "kc270"),
    ki270 = sum(pot == "ki270"),
    cc270 = sum(pot == "cc270"),
    ci270 = sum(pot == "ci270"),
    kc276 = sum(pot == "kc276"),
    ki276 = sum(pot == "ki276"),
    cc276 = sum(pot == "cc276"),
    ci276 = sum(pot == "ci276"),
    # co276 = sum(pot == "co276"),
  )

# add totals per row
surv_summary <- surv_summary %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("kc") | starts_with("ki") | starts_with("cc") | starts_with("ci") | starts_with("co")))) %>%
  ungroup()

# add totals per column
surv_summary <- surv_summary %>%
  adorn_totals("row")

# separate the last row (totals row)
last_row <- tail(surv_summary, 1)
surv_summary <- head(surv_summary, -1)

# apply conditional formatting to the pot count columns
surv_summary <- surv_summary %>%
  mutate(across(
    starts_with("kc") | starts_with("ki") | starts_with("cc") | starts_with("ci") | starts_with("co"),
    ~ cell_spec(.x, "html", color = "white", background = scales::col_numeric("Blues", domain = NULL)(.x))
  ))

# re-bind the last row back to the formatted data frame (with some type type tweaking...)
last_row <- last_row %>% mutate_all(as.character)
last_row$total <- as.integer(last_row$total)
surv_summary <- bind_rows(surv_summary, last_row)

# create and save pretty table
pretty_table <- kable(surv_summary, format = "html", escape = FALSE, table.attr = "class='table table-striped'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

pretty_table <- kable(surv_summary, format = "html", escape = FALSE, table.attr = "class='table table-striped'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  add_header_above(c("Infectious Bloodmeal" = 3, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, " " = 1)) %>%
  add_header_above(c("Mosquito species" = 3, "An. gambiae" = 2, "An. coluzzii" = 2, "An. gambiae" = 2, "An. coluzzii" = 2, "An. gambiae" = 2, "An. coluzzii" = 2, "An. gambiae" = 2, "An. coluzzii" = 2, " " = 1)) %>%
  add_header_above(c("Temperature range" = 3, "±0°C" = 4, "±6°C" = 4, "±0°C" = 4, "±6°C" = 4, " " = 1)) %>%
  add_header_above(c("Mean temperature" = 3, "21°C" = 8, "27°C" = 8, " " = 1))

pretty_table
save_kable(pretty_table, file = "Figures/sample_sizes.html")
webshot("Figures/sample_sizes.html", file = "Figures/sample_sizes.png", vwidth = 1400)

################################################################################################
# pretty table to show sample sizes and prevalence by pot and replicate
################################################################################################

# full disclosure - most of this stuff was written by a robot

pot_cols <- c(
  "kc210", "ki210", "cc210", "ci210",
  "kc216", "ki216", "cc216", "ci216",
  "kc270", "ki270", "cc270", "ci270",
  "kc276", "ki276", "cc276", "ci276" # , "co276"
)

# sample sizes
surv_summary_counts <- surv %>%
  group_by(replicate) %>%
  summarise(
    infection_date = first(date),
    length_days = as.character(max(dpi)),
    kc210 = sum(pot == "kc210"), ki210 = sum(pot == "ki210"),
    cc210 = sum(pot == "cc210"), ci210 = sum(pot == "ci210"),
    kc216 = sum(pot == "kc216"), ki216 = sum(pot == "ki216"),
    cc216 = sum(pot == "cc216"), ci216 = sum(pot == "ci216"),
    kc270 = sum(pot == "kc270"), ki270 = sum(pot == "ki270"),
    cc270 = sum(pot == "cc270"), ci270 = sum(pot == "ci270"),
    kc276 = sum(pot == "kc276"), ki276 = sum(pot == "ki276"),
    cc276 = sum(pot == "cc276"), ci276 = sum(pot == "ci276"),
    # co276 = sum(pot == "co276"),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(total = sum(c_across(all_of(pot_cols)))) %>%
  ungroup()

# prevalences
surv_summary_prev <- surv %>%
  group_by(replicate) %>%
  summarise(
    infection_date = first(date),
    length_days = as.character(max(dpi)),
    .groups = "drop"
  ) %>%
  left_join(
    surv %>%
      group_by(replicate, pot) %>%
      summarise(prev = round(mean(prevalence_oo, na.rm = TRUE), 3), .groups = "drop") %>%
      tidyr::pivot_wider(
        names_from = pot, values_from = prev,
        values_fill = 0
      ) %>%
      select(replicate, any_of(pot_cols)),
    by = "replicate"
  ) %>%
  # ensure all pot columns exist even if a pot never appears in any replicate
  tibble::add_column(!!!setNames(
    lapply(setdiff(pot_cols, names(.)), function(p) 0), # ← 0 not NA_real_
    setdiff(pot_cols, names(.))
  )) %>%
  select(replicate, infection_date, length_days, all_of(pot_cols)) %>%
  mutate(across(all_of(pot_cols), ~ replace(.x, is.na(.x) | is.nan(.x), 0))) %>% # ← catch any remaining
  rowwise() %>%
  mutate(mean = round(mean(c_across(all_of(pot_cols))[c_across(all_of(pot_cols)) != 0], na.rm = TRUE), 3)) %>%
  ungroup()
# view(surv_summary_prev)


# total midguts dissected (useful later)
surv_summary_midg <- surv %>%
  group_by(replicate, pot) %>%
  summarise(midgut_n = mean(midgut_n, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = pot, values_from = midgut_n, values_fill = 0) %>%
  select(replicate, any_of(pot_cols))

# view(surv_summary_midg)

# old code:
# surv_summary_prev <- surv_summary_prev %>%
#   mutate(
#     gametocytaemia = c(0.2, 0.2, 1.0, 0.3, 1.8, 0.8, 0.8, 1.2, 3.5, 2.7),
#     exflagellation = c("No", "No", "No", "No", "No", "No", "No", "No", "No", "No"),
#   )


surv_summary_prev <- surv_summary_prev %>%
  mutate(
    gametocytaemia = c(0.2, 0.2, 1.0, 0.3, 1.8, 0.8, 0.8, 1.2, 3.5, 2.7),
    age_gambiae = c(3.6, 3.6, 2.3, 1.8, 4.2, 3.1, 4.0, 5.0, 3.7, 3.5),
    age_coluzzii = c(2.9, NA, 2.2, 1.6, 3.5, 0.8, 3.3, 5.4, 3.7, 4.0)
  )

# function to build kable html string from dataframe
build_kable_html <- function(df, add_totals_row = TRUE, color_domain = NULL, totals_df = NULL, trow_name = NULL, is_prev = FALSE) {
  # old code
  # if (is_prev) {
  #   df <- df %>%
  #     select(-infection_date, -length_days) %>%
  #     rename(`gametocytaemia (%)` = gametocytaemia) %>%
  #     relocate(`gametocytaemia (%)`, exflagellation, .after = replicate)
  # }


  if (is_prev) {
    df <- df %>%
      select(-infection_date, -length_days) %>%
      rename(
        `G-taemia (%)` = gametocytaemia,
        `age (AG)` = age_gambiae,
        `age (AC)` = age_coluzzii
      ) %>%
      relocate(`G-taemia (%)`, `age (AG)`, `age (AC)`, .after = replicate)
  }

  body <- df

  if (add_totals_row) {
    last_row <- df %>%
      adorn_totals("row") %>%
      tail(1) %>%
      mutate_all(as.character)
    last_row$total <- as.numeric(last_row$total)
  }

  # if an external totals_df is supplied, build the bottom row from it
  if (!is.null(totals_df)) {
    if (is.null(trow_name)) {
      trow_name <- "Total"
    }

    ext_totals <- totals_df %>%
      adorn_totals("row", name = trow_name) %>%
      tail(1) %>%
      mutate_all(as.character)
    # keep only columns that exist in df, fill missing ones with ""
    ext_totals <- ext_totals %>%
      select(any_of(names(df))) %>%
      bind_cols(
        tibble::tibble(!!!setNames(
          lapply(setdiff(names(df), names(ext_totals)), function(x) "-"),
          setdiff(names(df), names(ext_totals))
        ))
      ) %>%
      select(all_of(names(df))) # restore column order
  }

  body_fmt <- body %>%
    mutate(
      # across(any_of(c("infection_date", "exflagellation")), as.character), # old code
      across(any_of("infection_date"), as.character),
      across(
        all_of(pot_cols),
        ~ cell_spec(
          replace(.x, is.na(.x) | is.nan(.x), 0),
          "html",
          color = "white",
          background = ifelse(
            replace(.x, is.na(.x) | is.nan(.x), 0) == 0,
            "white",
            scales::col_numeric("Blues", domain = color_domain)(
              replace(.x, is.na(.x) | is.nan(.x), 0)
            )
          )
        )
      )
    )

  if (add_totals_row) {
    body_fmt <- body_fmt %>% mutate(total = as.character(total))
    body_fmt <- bind_rows(body_fmt, last_row %>% mutate_all(as.character))
  }

  # append the external totals row (uncoloured, plain text) ──
  if (!is.null(totals_df)) {
    body_fmt <- body_fmt %>% mutate(across(where(is.numeric), as.character))
    body_fmt <- bind_rows(body_fmt, ext_totals %>% mutate_all(as.character))
  }

  kable(body_fmt,
    format = "html", escape = FALSE,
    table.attr = "class='table table-striped'"
  ) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    add_header_above(c(
      "Infectious bloodmeal" = ncol(df) - 17,
      "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1,
      "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, "No" = 1, "Yes" = 1, " " = 1
    )) %>%
    add_header_above(c(
      "Mosquito species" = ncol(df) - 17,
      "An. gambiae" = 2, "An. coluzzii" = 2, "An. gambiae" = 2, "An. coluzzii" = 2,
      "An. gambiae" = 2, "An. coluzzii" = 2, "An. gambiae" = 2, "An. coluzzii" = 2, " " = 1
    )) %>%
    add_header_above(c("Temperature range" = ncol(df) - 17, "±0°C" = 4, "±6°C" = 4, "±0°C" = 4, "±6°C" = 4, " " = 1)) %>%
    add_header_above(c("Mean temperature" = ncol(df) - 17, "21°C" = 8, "27°C" = 8, " " = 1))
}

# render both tables to html strings
html_counts <- as.character(build_kable_html(surv_summary_counts, add_totals_row = TRUE, color_domain = c(0, 222)))
html_prev <- as.character(build_kable_html(df = surv_summary_prev, add_totals_row = FALSE, color_domain = c(-0.4, 0.75), totals_df = surv_summary_midg, trow_name = "Midguts", is_prev = TRUE))

# build toggle page
page <- paste0('
<!DOCTYPE html><html lang="en"><head>
  <meta charset="UTF-8"><title>Sample Sizes</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css">
  <style>
    body { font-family: "Segoe UI", sans-serif; background:#f5f7fa; padding:2rem; }
    h2   { font-weight:700; color:#1a2b3c; }
    .toggle-bar { display:inline-flex; border:2px solid #1a6fa8; border-radius:8px;
                  overflow:hidden; margin-bottom:1.5rem; }
    .toggle-bar button { padding:.45rem 1.2rem; font-size:.9rem; font-weight:600;
                         border:none; background:white; color:#1a6fa8; cursor:pointer; }
    .toggle-bar button.active { background:#1a6fa8; color:white; }
    .table-wrapper { background:white; border-radius:10px;
                     box-shadow:0 2px 12px rgba(0,0,0,.08);
                     padding:1.2rem 1rem; overflow-x:auto; }
    #view-counts { display:block; }
    #view-prev   { display:none;  }
    .badge-mode  { font-size:.75rem; background:#e8f2fa; color:#1a6fa8;
                   padding:.2rem .6rem; border-radius:20px; margin-left:.5rem; }
  </style>
</head><body>
  <h2>Dataset summary <span id="mode-badge" class="badge-mode">Mosquitoes</span></h2>
  <div class="toggle-bar">
    <button id="btn-counts" class="active" onclick="show(\'counts\')">&#x23; Mosquitoes</button>
    <button id="btn-prev"                  onclick="show(\'prev\')"  >% Oocysts</button>
  </div>
  <div class="table-wrapper">
    <div id="view-counts">', html_counts, '</div>
    <div id="view-prev">', html_prev, '</div>
  </div>
  <script>
    function show(mode) {
      document.getElementById("view-counts").style.display = mode==="counts" ? "block":"none";
      document.getElementById("view-prev").style.display   = mode==="prev"   ? "block":"none";
      document.getElementById("btn-counts").classList.toggle("active", mode==="counts");
      document.getElementById("btn-prev").classList.toggle("active",   mode==="prev");
      document.getElementById("mode-badge").textContent = mode==="counts" ? "Mosquitoes":"Oocysts";
    }
  </script>
</body></html>')

writeLines(page, "Figures/sample_sizes_toggle.html")

page_counts_only <- paste0('
<!DOCTYPE html><html lang="en"><head>
  <meta charset="UTF-8">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css">
  <style>body{font-family:"Segoe UI",sans-serif;background:white;padding:1rem;}
         .table-wrapper{overflow-x:auto;}</style>
</head><body><div class="table-wrapper">', html_counts, "</div></body></html>")

writeLines(page_counts_only, "Figures/sample_sizes_counts_only.html")
webshot("Figures/sample_sizes_counts_only.html", "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/samplesizes.png", vwidth = 1400)

page_prev_only <- paste0('
<!DOCTYPE html><html lang="en"><head>
  <meta charset="UTF-8">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css">
  <style>body{font-family:"Segoe UI",sans-serif;background:white;padding:1rem;}
         .table-wrapper{overflow-x:auto;}</style>
</head><body><div class="table-wrapper">', html_prev, "</div></body></html>")

writeLines(page_prev_only, "Figures/sample_sizes_prev_only.html")
webshot("Figures/sample_sizes_prev_only.html", "Figures/sample_sizes_prev.png", vwidth = 1400)
webshot("Figures/sample_sizes_prev_only.html", "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/prevs.png", vwidth = 1400)

webshot("Figures/sample_sizes_toggle.html", "Figures/sample_sizes_toggle.png", vwidth = 1400)

# create new dataframe summarising survival outcome per replicate:pot
##############################################################################

# group surv by replicate:pot
surv_pot_summary <- surv %>%
  group_by(replicate, pot, pot_unique) %>% # creates one subdataframe for each unique combination of these variables (pot_unique is replicate+pot, but this ensures it is carried over)
  summarise( # creates a new dataframe with one row per group. The columns will be as described below:
    species = species[1], # species is the same for all rows
    temp_mean = temp_mean[1], # mean_temp is the same for all rows
    temp_range = temp_range[1], # temp_range is the same for all rows
    exposed = exposed[1], # exposure is the same for all rows
    treatment = paste0(temp_mean[1], "±", temp_range[1], "°C"), # will be useful for plotting
    .groups = "drop"
  )

# remove r00 because it was only monitored for 15 days, not 30
surv_pot_summary <- surv_pot_summary %>%
  filter(replicate != "r00")

# add new variable max_span
surv_pot_summary$max_span <- numeric(nrow(surv_pot_summary))
for (i in 1:nrow(surv_pot_summary)) {
  pot <- surv_pot_summary$pot_unique[i] # get pot name
  surv_pot_summary$max_span[i] <- max(surv[surv$pot_unique == pot, "dpi"]) # maximum lifespan in survival data
  .groups <- "drop"
}


# generate survival estimates based on pot-wise KM
##############################################################################

# initialise new vars (mean survival and CI limits)
surv_pot_summary$km_RMST_0toMax <- numeric(nrow(surv_pot_summary)) # integral of KM curve from t_0 to t_max
surv_pot_summary$km_RMST_0toMax_upper <- numeric(nrow(surv_pot_summary)) # integral of KM curve from t_0 to t_max + std.err
surv_pot_summary$km_RMST_0toMax_lower <- numeric(nrow(surv_pot_summary)) # integral of KM curve from t_0 to t_max - std.err
surv_pot_summary$n_total <- numeric(nrow(surv_pot_summary)) # number of mosquitoes in this pot
surv_pot_summary$n_deaths <- numeric(nrow(surv_pot_summary)) # number of deaths observed
surv_pot_summary$n_censored <- numeric(nrow(surv_pot_summary)) # number of mosquitoes censored


for (i in 1:length(unique(surv_pot_summary$pot_unique))) { # for each pot represented in surv_pot_summary
  pot <- surv_pot_summary$pot_unique[i] # get pot name
  print(pot)

  # fit a km curve for the current pot
  kmfit_pot <- survfit(event ~ 1, data = surv[surv$pot_unique == pot, ])

  # Extract full list of survival probabilities and times to calculate integrals
  surv_probs <- kmfit_pot$surv
  surv_times <- kmfit_pot$time
  surv_std_err <- kmfit_pot$std.err

  # # There's one case that has Inf std err... Not sure why, but that gets propagated so let's overwrite with NA
  # surv_std_err[is.infinite(surv_std_err)] <- NA

  # Since the KM is a step function, the integral is just a sum of rectangular areas
  # Since step = 1, all rectangles have width = 1. The sum of areas = the sum of survival probabilities

  # Integral t_0 to t_max
  surv_pot_summary$km_RMST_0toMax[i] <- sum(surv_probs)
  surv_pot_summary$km_RMST_0toMax_upper[i] <- sum(surv_probs + surv_std_err)
  surv_pot_summary$km_RMST_0toMax_lower[i] <- sum(surv_probs - surv_std_err)
  surv_pot_summary$n_total[i] <- kmfit_pot$n
  surv_pot_summary$n_deaths[i] <- sum(kmfit_pot$n.event)
  surv_pot_summary$n_censored[i] <- sum(kmfit_pot$n.censor)
}

# format variables correctly
surv_pot_summary$replicate <- as.factor(surv_pot_summary$replicate)
surv_pot_summary$pot <- as.factor(surv_pot_summary$pot)
surv_pot_summary$pot_unique <- as.factor(surv_pot_summary$pot_unique)
surv_pot_summary$species <- as.factor(surv_pot_summary$species)
surv_pot_summary$temp_mean <- as.factor(surv_pot_summary$temp_mean)
surv_pot_summary$temp_range <- as.factor(surv_pot_summary$temp_range)

# print and save dataset
# view(surv_pot_summary)
write.csv(surv_pot_summary, file = "Data/surv_pot_summary.csv", row.names = FALSE)

# create and save pretty table
pretty_table <- kable(surv_pot_summary, format = "html", escape = FALSE, table.attr = "class='table table-striped'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
save_kable(pretty_table, file = "Figures/surv_pot_summary_data.html")
save_kable(pretty_table, file = "Figures/surv_pot_summary_data.png")
