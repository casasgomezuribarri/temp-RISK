# this script defines custom functions

################################################################################################
# model comparison
################################################################################################
# for quick lrt between two models
calc_lrt <- function(big_model, small_model, name) {
  lrt_stat <- 2 * (logLik(big_model) - logLik(small_model))
  df_diff <- attr(logLik(big_model), "df") - attr(logLik(small_model), "df")
  p_val <- pchisq(as.numeric(lrt_stat), df = df_diff, lower.tail = FALSE)

  data.frame(Term = name, Chisq = round(as.numeric(lrt_stat), 4), Df = df_diff, p_value = round(p_val, 5))
}


compare_parametric_fits <- function(data # a lifetable
                                    , time_var # the variable with times
                                    , event_var # variable with status  (conventions from survival package)
                                    , dists = c("genf", "gengamma", "weibull", "gompertz", "gamma", "llogis", "lnorm", "exp") # distributions to fit, default is all of them
                                    , funcs = c("survival", "hazard", "cumhaz") # types of plots to produce, default is all of them
                                    , plot_title = "Parametric Survival Model Fits") {
  # useful for later
  pretty_dict <- c(
    genf     = "Gen. F",
    gengamma = "Gen. Gamma",
    weibull  = "Weibull",
    gompertz = "Gompertz",
    gamma    = "Gamma",
    llogis   = "Log-Logistic",
    lnorm    = "Log-Normal",
    exp      = "Exponential"
  )

  # fit all models
  fit_list <- list()
  for (dist in dists) {
    fit_list[[dist]] <- tryCatch(
      flexsurvreg(Surv(data[[time_var]], data[[event_var]]) ~ 1, data = data, dist = dist),
      error = function(e) {
        message(paste("Skipping", dist, "due to error:", e$message))
        NULL
      }
    )
  }

  # remove NULLS
  fit_list <- fit_list[!sapply(fit_list, is.null)]
  if (length(fit_list) == 0) { # dont continue if all NULL
    warning("No models were successfully fitted.")
    return(list(comparison = data.frame(), fits = NULL, plot = NULL))
  }

  # extract AIC, loglik, and pretty name for each fit
  aics <- sapply(fit_list, AIC)
  logliks <- sapply(fit_list, function(fit) fit$loglik)

  # generating a pvalue based on coxsnell residuals
  ks_pvals <- sapply(fit_list, function(fit) { # for each fitted model
    cs <- coxsnell_flexsurvreg(fit) # get the coxsnell residuals
    km_fit <- survfit(Surv(cs$est, cs$status) ~ 1) # fit a km curve to those residuals (*reason below)

    # data from km_fit
    obs_times <- km_fit$time
    emp_surv <- km_fit$surv

    # data from a theoretical population with exponential hazard
    theo_surv <- exp(-obs_times)

    # are they the same thing?
    ks_pval <- ks.test(emp_surv, theo_surv)$p.value
  })
  # *reason for that km fit
  # - CS residuals should resemble data from an Exp(1) without censoring
  # - fitting a KM to it effectively handles the censoring present
  # - check https://search.r-project.org/CRAN/refmans/flexsurv/html/coxsnell_flexsurvreg.html

  # use the dicitonary from above
  dist_names <- names(aics) # this is a named vector and we can use it
  dist_labels <- pretty_dict[dist_names]
  dist_labels[is.na(dist_labels)] <- dist_names[is.na(dist_labels)]

  # compile in a table, then sort it
  comparison <- data.frame(
    Distribution = dist_labels,
    AIC = aics,
    LogLik = logliks,
    KS_p = ks_pvals,
    row.names = NULL
  )
  comparison <- comparison[order(comparison$AIC), ]

  # reorder everything else (useful later)
  fit_list_ordered <- fit_list[comparison$Distribution %>% match(dist_labels)]
  # dist_ordered <- dist_names[order(aics)]
  # dist_labels_ordered <- dist_labels[order(aics)]


  # now let's produce plots in a table as well:
  n_dists <- length(fit_list_ordered)
  n_funcs <- length(funcs)
  op <- par(no.readonly = TRUE) # save old par settings
  par(mfrow = c(n_funcs, n_dists), mar = c(2, 2, 6, 1), oma = c(4, 4, 6, 2)) # margins and layout

  # # titles (top row)
  # for (label in dist_labels_ordered) {
  #   plot.new()
  #   title(main = label, cex.main = 1.2)
  # }

  # the actual plots
  for (func in funcs) {
    for (dist in dist_names) {
      fit <- fit_list[[dist]]
      if (func == funcs[1]) {
        plot(fit,
          type = func, xlab = "", ylab = func,
          main = paste0(pretty_dict[[dist]], "\nAIC: ", round(aics[dist], 2), "\np: ", round(ks_pvals[dist], 3))
        )
      } else {
        plot(fit, type = func, xlab = "", ylab = func, main = NULL)
      }
    }
  }

  # extras
  mtext("Time", side = 1, outer = TRUE, line = 2)
  mtext("", side = 2, outer = TRUE, line = 2)
  mtext(plot_title, side = 3, outer = TRUE, line = 4, cex = 1.5)

  plot <- recordPlot()

  par(op) # reset plotting parameters

  invisible(list(fits = fit_list, comparison = comparison, plot = plot))
}


# get empirical KM curves for custom groups (specifically tailored to my use case)
get_km_data <- function(surv_df, nd, grouping_col = NA) {
  # these labels are for plotting downstream
  nd <- nd %>%
    mutate(
      treatment = paste0(temp_mean, "±", temp_range, "°C")
    )

  # get the unique grouping combinations present in nd
  if (!is.na(grouping_col)) {
    groups <- nd %>%
      select(!!sym(grouping_col), species, temp_range, temp_mean, treatment) %>%
      distinct()
  } else {
    groups <- nd %>%
      select(species, temp_range, temp_mean, treatment) %>%
      distinct()
  }

  # fetch survival data
  km_list <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, ]

    # filter surv to rows matching this group
    if (!is.na(grouping_col)) {
      sub <- surv_df %>%
        filter(
          species == g$species,
          temp_range == g$temp_range,
          temp_mean == g$temp_mean,
          .data[[grouping_col]] == g[[grouping_col]]
        )
    } else {
      sub <- surv_df %>%
        filter(
          species == g$species,
          temp_range == g$temp_range,
          temp_mean == g$temp_mean
        )
    }


    # skip silently if no data
    if (nrow(sub) == 0) {
      return(NULL)
    }

    # fit km and format as tibble
    km <- survfit(Surv(age, dead) ~ 1, data = sub)
    if (!is.na(grouping_col)) {
      tibble(
        time = c(0, km$time),
        est = c(1, km$surv),
        lcl = c(1, km$lower),
        ucl = c(1, km$upper),
        species = g$species,
        treatment = g$treatment,
        !!sym(grouping_col) := g[[grouping_col]]
      )
    } else {
      tibble(
        time = c(0, km$time),
        est = c(1, km$surv),
        lcl = c(1, km$lower),
        ucl = c(1, km$upper),
        species = g$species,
        treatment = g$treatment
      )
    }
  })

  # apply to all groups and combine
  bind_rows(km_list)
}


# This function counts the amount of colons in a string.
# It is used to order formula terms so that interactions are evaluated first
count_colons <- function(x) {
  matches <- gregexpr(":", x)[[1]]
  if (matches[1] == -1) {
    return(0)
  } else {
    return(length(matches))
  }
}

# This function computes LRT statistics for all terms in a model
# only used for coxme and flexsurvreg mdoels
perform_lrt <- function(full_model, full_formula, data) {
  terms <- attr(terms(full_formula), "term.labels")

  # coxme and flexsurvreg models support frailty terms in the formula.
  # check if frailty term is in the formula, they need to be handled separately
  if (inherits(full_model, "coxme")) {
    has_replicate <- any(grepl("1 \\| replicate/pot", terms))
  } else if (inherits(full_model, "flexsurvreg")) {
    has_replicate <- any(grepl("cluster(replicate)", terms))
  }

  # remove if present (will add manually later)
  if (has_replicate) {
    if (inherits(full_model, "coxme")) {
      terms <- terms[terms != "1 | replicate/pot"]
    } else if (inherits(full_model, "flexsurvreg")) {
      terms <- terms[terms != "cluster(replicate)"]
    }
  }

  # initialise empty dataframe to store results
  results <- data.frame(Variable = character(), LRT_Statistic = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

  for (term in terms) {
    # create the reduced formula by removing one term at a time
    reduced_terms <- setdiff(terms, term)

    # add frailty term back if it was present in full_formula
    if (has_replicate) {
      if (inherits(full_model, "coxme")) {
        reduced_formula <- as.formula(paste("event ~", paste(reduced_terms, collapse = " + "), "+ (1|replicate/pot)"))
      } else if (inherits(full_model, "flexsurvreg")) { # this is probably not necessary but oh well
        reduced_formula <- as.formula(paste("event ~", paste(reduced_terms, collapse = " + "), "+ cluster(replicate)"))
      }
    } else {
      reduced_formula <- as.formula(paste("event ~", paste(reduced_terms, collapse = " + ")))
    }

    # before fitting model, check if it will return an error (usually bc non-convergence)
    reduced_model <- tryCatch(
      {
        if (inherits(full_model, "coxme")) {
          coxme(reduced_formula, data = data)
        } else if (inherits(full_model, "flexsurvreg")) {
          flexsurvreg(reduced_formula, data = data, dist = full_model$dlist$name)
        } else {
          stop("Unsupported model type (must be coxme or flexsurvreg)")
        }
      },
      error = function(e) {
        return(list(error = e$message))
      }
    )

    # if it doesn't converge, store error message in lrt column so can be checked. Then p <- NA
    if (is.list(reduced_model) && !is.null(reduced_model$error)) {
      p_value <- NA
      lrt_statistic <- reduced_model$error
    } else {
      # Calculate the LRT statistic and p-value
      logLik_full <- logLik(full_model)
      logLik_reduced <- logLik(reduced_model)
      lrt_statistic <- -2 * (logLik_reduced - logLik_full)
      df <- abs(attr(logLik_full, "df") - attr(logLik_reduced, "df")) # is this okay? sometimes df is marginally negative (suuuuper low |df|...)
      p_value <- pchisq(lrt_statistic, df, lower.tail = FALSE)
    }
  }

  # Store the results
  results <- rbind(results, data.frame(Variable = term, LRT_Statistic = lrt_statistic, p_value = p_value))
  return(results)
}

# this function removes level from variable name and deals with the order of dpi when interacting (needs attention though)
clean_term_name <- function(term) {
  term <- gsub("exposedExposed|exposedControl", "exposed", term)
  term <- gsub("speciesAn. gambiae|speciesAn. coluzzii", "species", term)
  term <- gsub("infectedUninfected|infectedInfected", "infected", term)
  term <- gsub("temp_mean27|temp_mean21", "temp_mean", term)
  term <- gsub("treatment27 ± 0|treatment21 ± 0|treatment27 ± 6|treatment21 ± 6", "treatment", term)
  term <- gsub("temp_range0|temp_range6", "temp_range", term)
  term <- gsub("treatment21 ± 6|treatment21 ± 0|treatment27 ± 6|treatment27 ± 0", "treatment", term)
  term <- gsub("dpi:(.*)", "\\1:dpi", term) # if dpi appears at the beginning, make it appear at the end of the interaction instead - honestly why does it changeee
  return(term)
}
# this is super dirty and too specfic. This needs to be updated so that I can extract the VARIABLE WITHOUT THE LEVEL
# must work with interactions: "temp_mean27:speciesAn. gambiae"  |> "temp_mean:species"
# why does the order change when interacting with dpi?

# this function extracts the relevant info from a coxphw model and presents in a better way - we might not need aaaaall of this
extract_coxphw_results <- function(model, formula) {
  terms <- attr(terms(formula), "term.labels") # list of variables in the model
  coefs <- model$coefficients # this is similar, but includes the relevant level if a factor
  results <- data.frame(
    var = character() # intialise empty dataframe
    , clean_var = character(),
    coef = numeric(),
    hr = numeric(),
    hr_lower = numeric(),
    hr_upper = numeric(),
    p_coef = numeric(),
    Wald_Statistic = numeric(),
    p_wald = numeric()
  )

  for (var in names(coefs)) {
    coef_var <- coefs[var] # predictor of this coefficient
    hr_var <- exp(coef_var)
    ci_low <- model$ci.lower[var]
    ci_hig <- model$ci.upper[var]
    p <- model$prob[var]
    var_varcov <- model$var[var, var, drop = FALSE]
    wald_result <- wald(coef_var, var_varcov)
    results <- rbind(
      results, data.frame(
        var = var # predictor, but without the factor level
        , clean_var = clean_term_name(var),
        coef = coef_var # coefficient
        , hr = hr_var # HR
        , hr_lower = ci_low # upper ci of hr
        , hr_upper = ci_hig # lower ci of hr
        , p_coef = p # p for coefficient
        , Wald_Statistic = wald_result["chi2"] # wald statistic
        , p_wald = wald_result["p"]
      ) # p for wald statistic
    )
  }

  return(results)
}

# handle random effects in coxme and flexsurvreg. I don't think this functino does what is supposed to... will come back to this later.
create_reduced_formula <- function(terms, model) {
  if (inherits(model, "coxme")) {
    has_replicate <- any(grepl("1 \\| replicate", terms)) # checks if it has replicate as a random effect
    if (has_replicate) {
      terms <- terms[terms != "1 | replicate"] # removes it
      return(as.formula(paste("event ~", paste(terms, collapse = " + "), "+ (1|replicate)"))) # adds it in correct format
    }
  } else if (inherits(model, "flexsurvreg")) { # same as above
    has_replicate <- any(grepl("cluster(replicate)", terms))
    if (has_replicate) {
      terms <- terms[terms != "cluster(replicate)"]
      return(as.formula(paste("event ~", paste(terms, collapse = " + "), "+ cluster(replicate)")))
    }
  }
  return(as.formula(paste("Surv(dpi, dead) ~", paste(terms, collapse = " + "))))
}

# a tryctch block
try_fit_model <- function(formula, model, data) {
  tryCatch(
    {
      if (inherits(model, c("coxme", "flexsurvreg", "coxphw"))) {
        new_call <- model$call
        new_call$formula <- formula
        return(eval(new_call)) # fit model with all the same arguments but the new formula
      } else {
        stop("Unsupported model type. Must be coxme, flexsurvreg or coxphw")
      }
    },
    error = function(e) {
      return(list(error = e$message))
    }
  )
}

# helper function to deal with interaction terms that are unordered
are_terms_equivalent <- function(term1, term2) {
  split_term1 <- sort(unlist(strsplit(term1, ":")))
  split_term2 <- sort(unlist(strsplit(term2, ":")))
  return(identical(split_term1, split_term2))
}

# augmented backward elimination
find_best_model_abe <- function(full_model, full_formula, data, active = NULL, passive = NULL, p_thr = 0.2, t_thr = 0.05, logfile = "model_log.txt") { # "active" and "passive" variable labels default to NULL but can be parsed.
  if (!file.exists(logfile)) {
    file.create(logfile)
  }
  sink(logfile)

  print(Sys.time())
  print("Curent test: ")
  print(" - formula:")
  print(full_formula)
  print(paste(" - model:", full_model$call))
  print(paste(" - method: augmented backwards elimination"))
  print(paste(" - p:", p_thr))
  print(paste(" - tau:", t_thr))
  if (exists("fun")) { # this is not robust... But works for now.
    fun_definition <- paste(deparse(fun), collapse = "\n")
    print(paste(" - f(t):", fun_definition))
  } else {
    # Log an alternative message
    print(" - f(t): t")
  }
  print("           ")


  best_model <- full_model
  best_formula <- full_formula

  repeat {
    # extract the model's output
    if (inherits(best_model, "coxphw")) {
      results <- extract_coxphw_results(best_model, best_formula)
      active_and_passive <- setdiff(attr(terms(best_formula), "term.labels"), c(active, passive)) # all variables are active_and_passive unless otherwise specified
    } else {
      results <- perform_lrt(best_model, full_formula, data) # needs to be updated
    }
    # sort results according to p_coef and colon number (interactions first, then high p-values first)
    results <- results[order(-sapply(results$var, count_colons), -results$p_coef), ]

    # evaluation of significance: if p > p_thr AND var is not "passive", consider for elimination
    black_list <- results[!results[, "clean_var"] %in% passive & results$p_coef > p_thr, ]
    rownames(black_list) <- NULL # important..

    # print msg:
    print("Active variables based on significance:")
    print(black_list[, c("var", "p_coef")])
    print("        ")

    # if black_list is empty, break the loop
    if (nrow(black_list) == 0) {
      print("No more variables to remove based on significance.")
      break
    }

    # evaluation of change-in-estimate criterion
    term_removed <- FALSE
    for (act_var in black_list[, "var"]) { # iterator through variables in black_list
      # act_var = black_list[, "var"][1]
      act_var_idx <- which(results[, "var"] == act_var)
      term_to_remove <- results[, "clean_var"][act_var_idx] # it's useful to have it in this form too
      print(paste("Evaluating", act_var))

      # if term_to_remove exists in best_formula interacting with dpi, we can't remove it!
      interaction_exists <- FALSE # default
      interaction_terms <- c(paste0(term_to_remove, ":dpi"), paste0(term_to_remove, ":fun(dpi)"), paste0("dpi:", term_to_remove), paste0("fun(dpi):", term_to_remove))
      for (int_term in interaction_terms) {
        if (int_term %in% attr(terms(best_formula), "term.labels")) {
          print("Found this term interacting with time, can't remove it!")
          print("             ")
          interaction_exists <- TRUE
          break
        }
      }
      if (interaction_exists) {
        next # move to next act_var
      }

      coef_act <- results[act_var, "coef"] # coefficient of active variable
      var_aa <- best_model$var[act_var, act_var, drop = FALSE] # variance of active coefficient
      skip_to_next_act_var <- FALSE # useful later
      cies <- data.frame(
        active = character() # initialise a dataframe for logs
        , passive = character(),
        scaled_cie = numeric(),
        p_active = numeric()
      )

      for (pas_var_temp in intersect(c(active_and_passive, passive), results[, "clean_var"])) { # iterator through variables with "passive" label that are present in the current model
        # pas_var_temp = intersect(c(active_and_passive, passive), results[, "clean_var"])[1]
        term_to_compare <- pas_var_temp
        pas_var_indices <- which(results[, "clean_var"] == pas_var_temp) # (could be more than 1 if pas_var is a factor with >2 levels)
        for (pas_var_idx in pas_var_indices) {
          # pas_var_idx = pas_var_indices[1]
          pas_var <- results[pas_var_idx, "var"] # overwrite with ugly name (dummy variable)

          if (pas_var == act_var | term_to_compare == term_to_remove) { # don't compare it to itself... Obviously removing a variabel changes its coefficient!
            next
          }

          var_pa <- best_model$var[act_var, pas_var, drop = FALSE] # covariance of passive and active estimates

          # if pas_var is an interaction with time, sd is not reported and needs to be calculated manually
          # we restrict a bit the way in which interactions with time can be specified... see cases below and error message.
          if (grepl("fun(dpi):", pas_var)) {
            rest <- gsub("fun\\(dpi\\):", "", pas_var) # get the term interacting with time
            sd_pp <- sd(model.matrix(best_model)[, rest] * data$dpi) # manually calculate the sd
          } else if (grepl("dpi:", pas_var)) {
            rest <- gsub("dpi:", "", pas_var) # get the term interacting with time
            sd_pp <- sd(model.matrix(best_model)[, rest] * data$dpi) # manually calculate the sd
          } else if (grepl(":fun(dpi)", pas_var)) {
            rest <- gsub(":fun\\(dpi\\)", "", pas_var) # get the term interacting with time
            sd_pp <- sd(model.matrix(best_model)[, rest] * data$dpi) # manually calculate the sd
          } else if (grepl(":dpi", pas_var)) {
            rest <- gsub(":dpi", "", pas_var) # get the term interacting with time
            sd_pp <- sd(model.matrix(best_model)[, rest] * data$dpi) # manually calculate the sd
          } else if (grepl("dpi", pas_var)) {
            stop("Interaction with time must be specified as 'var:dpi' or 'var:fun(dpi)'. Functions to it need to be defined elsewhere and should be called fun().")
          } else {
            sd_pp <- sd(model.matrix(best_model)[, pas_var]) # standard deviation of passive variable (model.matrix ensures that if pas_var is a factor, the appropriate dummy var is considered)
          }
          # to do: add function argument specifying time variable - won't always be dpi! - this'll work for me though.
          # or fina way to extract its name from the function call! - this may differ across model classes though. here's a starting point:
          # model_formula <- formula(full_model)
          # time_var <- all.vars(model_formula[[2]])[1]
          # print(time_var)

          change_in_estimate <- -(coef_act * var_pa / var_aa) # as in Dunkler et al 2014
          cie_scaled <- abs(change_in_estimate) * sd_pp # as in Dunkler et al 2014
          if (cie_scaled >= log(1 + t_thr)) { # this is the cie criterion test. If true, cie too big, active var stays.
            print("               ")
            print("Change-in-estimate criterion not satisfied: ")
            print(paste("Removing", act_var, "causes substantial change in coefficient of", pas_var, "cie =", cie_scaled))
            print("                                                  ")
            skip_to_next_act_var <- TRUE # create variable to break second outer loop
            break # if any passive variable changes more than t_thr when removing current active one, go to next in blacklist, current active var stays in the model
          }

          # log results
          new_row <- data.frame(
            active = act_var # active variable
            , passive = pas_var # passive variable
            , scaled_cie = cie_scaled[1] # change in estimate, scaled by sd(x_p)
          )
          rownames(new_row) <- NULL # important..


          cies <- rbind(cies, new_row)
          rownames(cies) <- NULL # important..

          print(paste("Change in", pas_var, ":", cie_scaled[1]))
        }
        if (skip_to_next_act_var) { # if we exited the last loop because we failed the cie test, break loop over passive variables
          break
        }
      }
      if (skip_to_next_act_var) { # if we exited the last loop because we failed the cie test, break loop and jump to next active variable
        next
      }

      # remove term from formula

      current_formula_terms <- attr(terms(best_formula), "term.labels")
      matching_term <- list() # there should always be one and only one

      for (var in current_formula_terms) {
        if (are_terms_equivalent(var, act_var)) {
          matching_term <- append(matching_term, list(var))
        }
      }

      reduced_terms <- setdiff(current_formula_terms, matching_term)
      reduced_formula <- create_reduced_formula(reduced_terms, full_model) # full_model here is used to check model type

      # fit new model
      print(paste("Removing term", term_to_remove, "and fitting new model with formula:"))
      print(paste("                                                  "))
      print(reduced_formula)
      best_model <- try_fit_model(reduced_formula, full_model, data) # copy model call from full_model but update formula to reduced_formula
      if (is.null(best_model$error)) {
        best_formula <- reduced_formula
        term_removed <- TRUE
        break
      } else {
        print(paste("Couldn't fit the following model:", best_model))
        print(best_model$error)
        break
      }
    }

    if (!term_removed) {
      print("No terms were removed in this iteration")
      break
    }
  }

  sink() # break connection with log file
  return(list(best_model = best_model, best_formula = best_formula))
}


################################################################################################
# plotting
################################################################################################

# km  curves

create_survival_curve <- function(data, title, colors, title_size = 35, group = "exposed", annotate = TRUE) {
  if (!(is.na(group))) {
    # convert group to a factor, define formula, fit KM model
    data[[group]] <- as.factor(data[[group]])
    surv_formula <- as.formula(paste("data$event ~ data$", group, sep = ""))
    fit <- survfit2(surv_formula)

    # extract median survival for each level in data[[group]]
    surv_data <- summary(fit)$table
    annotations <- data.frame(
      group = rownames(surv_data),
      x = surv_data[, "median"]
    )
    labels <- levels(data[[group]]) # labels for the annotation


    # plot
    p <- survfit2(surv_formula) %>%
      ggsurvfit() +
      labs(
        x = "Time (days)",
        y = "Proportion alive"
      ) +
      xlim(0, 35) + # Make sure all plots have the same axis limits
      ylim(0, 1) + # Make sure all plots have the same axis limits
      add_confidence_interval() +
      scale_fill_manual(values = colors) + # Set fill colors based on the palette
      scale_color_manual(values = colors) + # Set line colors based on the palette
      theme(
        axis.text = element_text(size = 30), # Font size for axis ticks
        axis.title = element_text(size = 40), # Adjust font of labels
        legend.position = "none", # Try "bottom"
        plot.margin = margin(10, 10, 10, 10), # Plot margins (t, r, b, l)
        plot.title = element_text(size = title_size, hjust = 0.5) # Title settings
      ) +
      ggtitle(title)

    if (annotate == TRUE) {
      # add a single "Median survival" annotation
      p <- p + annotate("text",
        x = 0, y = 0.25,
        label = "Median survival:",
        size = 8, color = "black", fontface = "bold", hjust = 0
      )

      # add the median survival for each level, and the lines
      for (i in 1:nrow(annotations)) {
        # add horizontal line
        p <- p + geom_segment(
          x = 0, xend = annotations$x[i], y = 0.5, yend = 0.5,
          color = colors[i], linetype = "solid"
        )

        # add vertical line
        p <- p + geom_segment(
          x = annotations$x[i], xend = annotations$x[i], y = 0.5, yend = 0,
          color = colors[i], linetype = "solid"
        )

        # Add annotation text
        p <- p + annotate("text",
          x = 0, y = 0.25 - i * 0.05, # Adjust y position for each label
          label = paste(labels[i], ": ", round(annotations$x[i], 1)),
          size = 8, color = colors[i], hjust = 0
        )
      }
    }
  } else {
    # convert group to a factor, define formula, fit KM model
    surv_formula <- as.formula("data$event ~ 1")
    fit <- survfit2(surv_formula)

    # extract median survival for each level in data[[group]]
    surv_data <- summary(fit)$table
    annotations <- data.frame(
      group = c("All"),
      x = surv_data["median"]
    )
    labels <- c("All")


    # plot
    p <- survfit2(surv_formula) %>%
      ggsurvfit() +
      labs(
        x = "Time (days)",
        y = "Proportion alive"
      ) +
      xlim(0, 35) + # Make sure all plots have the same axis limits
      ylim(0, 1) + # Make sure all plots have the same axis limits
      add_confidence_interval() +
      scale_fill_manual(values = colors) + # Set fill colors based on the palette
      scale_color_manual(values = colors) + # Set line colors based on the palette
      theme(
        axis.text = element_text(size = 30), # Font size for axis ticks
        axis.title = element_text(size = 40), # Adjust font of labels
        legend.position = "none", # Try "bottom"
        plot.margin = margin(10, 10, 10, 10), # Plot margins (t, r, b, l)
        plot.title = element_text(size = title_size, hjust = 0.5) # Title settings
      ) +
      ggtitle(title)
  }
  return(p)
}
create_legend_plot <- function(
  data, colors, ltypes,
  group = "exposed"
) {
  # convert the group to a factor and define formula
  data[[group]] <- as.factor(data[[group]])
  surv_formula <- as.formula(paste("data$event ~ data$", group, sep = ""))

  # define labels
  if (group == "exposed") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Exposure:"
  } else if (group == "species") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Species:"
  } else if (group == "temp_mean") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Mean temperature:"
  } else if (group == "temp_range") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Temperature oscillations:"
  } else {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Group:"
  }
  # maybe double check in which order each should be... but it seems fine to me like this

  # first make a plot
  p <- survfit2(surv_formula) %>%
    ggsurvfit() +
    scale_fill_manual(values = colors, labels = labels) + # Set fill colors and labels
    scale_color_manual(values = colors, labels = labels) + # Set line colors and labels
    add_confidence_interval() +
    theme(
      legend.text = element_text(size = 35), # Adjust legend text size
      legend.position = "bottom", # Position the legend at the bottom
      legend.title = element_text(size = 35), # Adjust legend title size
    ) +
    guides(
      fill = guide_legend(title = l_title),
      color = guide_legend(title = l_title)
    )

  # then xtract the legend
  grob <- ggplotGrob(p) # convert the plot to a grob
  # print(grob$layout$name) # find where the legend is (guide-box-bottom)
  positions <- which(grob$layout$name == "guide-box-bottom")
  legend <- grob$grobs[[positions[1]]] # extract the  legend
  # as_ggplot(legend) # plot the legend
  return(legend)
}


# model outputs

plot_output <- function(model, title = "Hazard Ratios with 95% CI") {
  # extract useful data from model
  results <- data.frame(
    Variable = names(model$coefficients),
    HR = exp(model$coefficients),
    CI_Lower = model$ci.lower,
    CI_Upper = model$ci.upper,
    P_Value = model$prob
  )

  # plot
  p <- ggplot(results, aes(x = Variable, y = HR)) +
    geom_point() +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.08) +
    geom_text(aes(label = ifelse(P_Value < 0.001, "p<0.001", sprintf("p = %.3f", P_Value))), vjust = -1) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::label_number(accuracy = 0.001)
    ) +
    labs(
      title = title,
      y = "Hazard Ratio (log scale)",
      x = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 25, face = "bold"),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20)
    ) +
    coord_flip() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") # add a dashed line at x = 1

  print(p)
}


################################################################################################
# conditional survival
################################################################################################
# data = surv
# column = "pot"
# value = "kc270"
# time_var = "dpi"
# event_var = "dead"
# taus = seq(1, 30, 1) # seq(1, 30, 1) || c(5, 10, 20, 30)

# conf_int_level = 0.95


# conditional survival function - P(T>tau | T>t) ∀ tau ∈ taus, ∀ t ∈ [1, max_time]
compute_conditional_survival <- function(data, column, value, time_var, event_var, taus, conf_int_level = 0.95) {
  data_subset <- data[data[[column]] == value, ] # subset data based on column value

  # since this metric will be a ratio, we need to make sure that the propagated total uncertainty corresponds to the 95%CI
  ci <- 1 - 2 * sqrt((1 - conf_int_level) / 2) # setting ci to this value ensures that the final ratio will have 95%CI

  # first fit km curve for our subset
  kmfit <- survfit(Surv(data_subset[[time_var]], data_subset[[event_var]]) ~ 1,
    data = data_subset,
    conf.int = ci
  )

  # get the maximum age in the subset - conditional survival won't be defined beyond this.
  max_time <- max(data_subset[[time_var]], na.rm = TRUE)

  # define relevant time range
  time_points <- seq_len(max_time) # ℕ ∈ [1, max_time] - all the possible ages under current condition.
  A_t <- approx(x = kmfit$time, y = kmfit$surv, xout = time_points, yleft = 1)$y # Age pyramid, A(t) - survival function, for ℕ ∈ [1, max_time]
  A_t <- A_t / sum(A_t) # scaled to sum 1 (PMF)

  # compute conditional survival - here's the 'hing!
  results <- map_dfr(taus, function(tau) { # mapping > for loop (vectorised ops > rowwise ops). We'll do the following ∀ tau ∈ taus
    # tau = taus[1]

    # extract empirical P(T>tau), and its { 1-2*sqrt((1-conf_int_level)/2) } % CI limits from the km curve
    S_tau <- summary(kmfit, times = tau)$surv
    S_tau_lo <- summary(kmfit, times = tau)$lower
    S_tau_up <- summary(kmfit, times = tau)$upper

    # extract empirical P(T>t) ∀ t ∈ [1, max_time], and their { 1-2*sqrt((1-conf_int_level)/2) } % CI limits - also from the km curve
    S_t <- summary(kmfit, times = time_points)$surv
    S_t_lo <- summary(kmfit, times = time_points)$lower
    S_t_up <- summary(kmfit, times = time_points)$upper


    # initialise vectors to store results for this tau (default value is undefined)
    Stau_St <- rep(NA, length(time_points))
    Stau_St_lo <- rep(NA, length(time_points))
    Stau_St_up <- rep(NA, length(time_points))

    # we'll only populate results for valid indices (S(t) > 0 and tau <= max_time)
    valid_idx <- which(S_t > 0 & tau <= max_time)

    # compute P(T > tau | T > t) only for valid indices
    Stau_St[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau / S_t[valid_idx])
    Stau_St_lo[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau_lo / S_t_up[valid_idx])
    Stau_St_up[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau_up / S_t_up[valid_idx])

    # compute sum(S(t)) for all t > tau
    sum_St <- sum(S_t)
    sum_S_tau_end <- sum(S_t[time_points >= tau])

    # return a dataframe with results for this tau
    tibble(
      time = time_points,
      age_dist = A_t,
      tau = tau,
      Stau = S_tau, # this is useful for a little experiment
      sum_St = sum_St, # this is useful for a little experiment
      sum_S_tau_end = sum_S_tau_end, # this is useful for a little experiment
      Stau_St = Stau_St,
      Stau_St_lo = Stau_St_lo,
      Stau_St_up = Stau_St_up
    )
  })

  return(results)
}
