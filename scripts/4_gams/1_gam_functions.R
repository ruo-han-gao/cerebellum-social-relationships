# Function for running GAMs


###############################################
#####               Fit GAMs             ######
###############################################

fit.models <- function(formula, labels, dataset = dat) {
  model_list <- lapply(labels, function(region) {
    modelformula <- as.formula(sprintf("%s ~ %s", region, formula))
    mgcv::gam(modelformula, method = "REML", data = dataset)
  })
  names(model_list) <- labels
  return(model_list)
}


###############################################
#####             Compare GAMs           ######
###############################################

compare.models <- function(modlist1, modlist2) {
  
  # store results in a list, then bind later
  results <- vector("list", length(modlist1))
  
  # loop through each model
  for (i in seq_along(modlist1)) {
    mod1 <- modlist1[[i]]
    mod2 <- modlist2[[i]]
    
    # extract parcel name 
    dv.label <- all.vars(formula(mod1))[[1]]
  
    # get model fit indices
    aic1 <- AIC(mod1)
    aic2 <- AIC(mod2)
    
    # get r2
    mod1.r2 <- summary(mod1)$r.sq  
    mod2.r2 <- summary(mod2)$r.sq  
    
    anov <- anova(mod1, mod2, test = "F")
    F.value <- anov$F[2]
    P.value <- anov$`Pr(>F)`[2]
    
    results[[i]] <- data.frame(
      DV = dv.label,
      AIC.dif = round(aic1 - aic2, 3),
      Mod1_R2 = round(mod1.r2, 3),
      Mod2_R2 = round(mod2.r2, 3),
      delta_R2 = round(mod2.r2 - mod1.r2, 3),
      F.value = round(F.value, 3),
      P.value = round(P.value, 3))
  }
  
  model_compare <- do.call(rbind, results)
  rownames(model_compare) <- NULL
  return(model_compare)
}

better.models <- function(model.compare.df) {
  # Identify parcels where the extended model is better
  better.parcels <- model.compare.df %>%
    filter((P.value < 0.05 & AIC.dif > 1) | AIC.dif > 5) %>%
    pull(DV)
  
  # Calculate and print percentage of better parcels
  perc <- length(better.parcels) / nrow(model.compare.df) * 100
  cat("Percentage of better parcels:", round(perc, 1), "%\n")
  cat("Better parcels:", paste(better.parcels, collapse = ", "), "\n")
  
  return(better.parcels)
}

###############################################
#####      Fit GAMs for social vars      ######
###############################################

fit_hirac_models <- function(predictor, baseline_model, labels) {
  
  # Prefix for naming
  if (startsWith(predictor, "raw_")) {
    prefix_full <- substr(predictor, 5, nchar(predictor))  # remove "raw_"
  } else {
    prefix_full <- predictor
  }
  prefix <- substr(prefix_full, 1, 3)  # keep max 3 letters
  
  # Smooth-only model
  mod_name <- paste0("mod.", prefix, "smoo")
  mod_smoo <- fit.models(
    formula = paste0("s(", predictor, ") + s(icv) + s(age) + sex"),
    labels = labels
  )
  compare_name <- paste0("compare.", prefix, "smoo")
  compare_smoo <- compare.models(baseline_model, mod_smoo)
  best_name <- paste0(prefix, "smoo.better")
  smoo_better <- better.models(compare_smoo)
  
  # Linear age interaction
  mod_age1_name <- paste0("mod.", prefix, "age1")
  mod_age1 <- fit.models(
    formula = paste0(predictor, ":age + s(", predictor, ") + s(icv) + s(age) + sex"),
    labels = labels
  )
  compare_age1_name <- paste0("compare.", prefix, "age1")
  compare_age1 <- compare.models(baseline_model, mod_age1)
  best_age1_name <- paste0(prefix, "age1.better")
  age1_better <- better.models(compare_age1)
  
  # Tensor age interaction (nonlinear)
  mod_age2_name <- paste0("mod.", prefix, "age2")
  mod_age2 <- fit.models(
    formula = paste0("ti(", predictor, ", age) + s(", predictor, ") + s(icv) + s(age) + sex"),
    labels = labels
  )
  compare_age2_name <- paste0("compare.", prefix, "age2")
  compare_age2 <- compare.models(mod_age1, mod_age2)
  best_age2_name <- paste0(prefix, "age.better")
  age2_better <- better.models(compare_age2)
  
  # Sex interaction
  mod_sex_name <- paste0("mod.", prefix, "sex")
  mod_sex <- fit.models(
    formula = paste0("s(", predictor, ", by = sex) + s(icv) + s(age) + sex"),
    labels = labels
  )
  compare_sex_name <- paste0("compare.", prefix, "sex")
  compare_sex <- compare.models(baseline_model, mod_sex)
  best_sex_name <- paste0(prefix, "sex.better")
  sex_better <- better.models(compare_sex)
  
  # Return all with dynamic names
  return(list(
    models = setNames(list(mod_smoo, mod_age1, mod_age2, mod_sex),
                      c(mod_name, mod_age1_name, mod_age2_name, mod_sex_name)),
    comparisons = setNames(list(compare_smoo, compare_age1, compare_age2, compare_sex),
                           c(compare_name, compare_age1_name, compare_age2_name, compare_sex_name)),
    best = setNames(list(smoo_better, age1_better, age2_better, sex_better),
                    c(best_name, best_age1_name, best_age2_name, best_sex_name))
  ))
}


###############################################
#####           stats from GAMs          ######
###############################################

get.gam.stats <- function(modelist, model_name = deparse(substitute(modelist))) {
  
  mod.stat <- lapply(modelist, function(mod) {
    
    # get model summary
    summ <- summary(mod)
    
    process_table <- function(df) {
      df <- as.data.frame(df)
      df$term <- rownames(df)
      
      # name the columns as gam.stat/lm.stat
      names(df)[1:4] <- c("edf/estimate", "Ref.df/se", "F/t", "p.value")
      df
    }
    # get stats for smooth terms
    gam.stat <- process_table(summ$s.table)
    
    # get stats for linear terms
    lm.stat  <- process_table(summ$p.table)
    
    # combine them
    stat <- bind_rows(gam.stat, lm.stat)
    # add the passed model list name
    stat$model_name <- model_name
    
    # Add region and R2
    stat$DV <- all.vars(formula(mod))[[1]]
    stat$r.squr <- summ$r.sq
    
    return(stat)
  })
  
  # combine all 32 models
  all.stat <- bind_rows(mod.stat)
  terms_out_interest <- c("s(icv)", "s(age)", "(Intercept)", "sexF_M")
  
  # apply FDR across all models per term
  all.stat <- all.stat %>%
    group_by(term) %>%
    mutate(
      p.fdr = p.adjust(p.value, method = "fdr")
    ) %>%
    ungroup() %>%  # ungroup after computing p.fdr
    mutate(
      sig = case_when(
        !term %in% terms_out_interest & p.fdr < 0.05   ~ "YES",
        !term %in% terms_out_interest & p.value < 0.05 ~ "yes",
        term %in% terms_out_interest & p.fdr < 0.05   ~ "Y",
        term %in% terms_out_interest & p.value < 0.05 ~ "y",
        TRUE ~ NA_character_
      )
    )
  
  # Round numeric columns
  num_cols <- intersect(names(all.stat), c("edf/estimate", "Ref.df/se", "F/t", "p.value", "r.squr", "p.fdr"))
  all.stat[num_cols] <- lapply(all.stat[num_cols], function(x) round(as.numeric(x), 3))
  
  return(all.stat)
}

get.int.stats <- function(modelist, model_name = deparse(substitute(modelist))) {
  
  mod.stat <- lapply(modelist, function(mod) {
    
    summ <- summary(mod)
    
    process_table <- function(df) {
      df <- as.data.frame(df)
      df$term <- rownames(df)
      names(df)[1:4] <- c("beta_coef", "se", "t", "p.value")
      df
    }
    
    lm.stat <- process_table(summ$p.table)
    lm.stat$DV     <- all.vars(formula(mod))[[1]]
    
    # add the passed model list name
    lm.stat$model_name <- model_name
    
    return(lm.stat)
  })
  
  # combine all model outputs
  all.stat <- bind_rows(mod.stat)
  
  # remove intercept
  all.stat <- all.stat %>% filter(!term %in% c("(Intercept)", "sexM"))
  
  # compute FDR per term
  all.stat <- all.stat %>%
    group_by(term) %>%
    mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>%
    ungroup()
  
  # significance labels
  all.stat <- all.stat %>%
    mutate(
      sig = case_when(
        p.fdr < 0.05 ~ "YES",
        p.value < 0.05 ~ "yes",
        TRUE ~ NA_character_
      )
    )
  
  # round numeric columns
  num_cols <- c("beta_coef", "se", "t", "p.value", "p.fdr")
  all.stat[num_cols] <- lapply(all.stat[num_cols], function(x) round(as.numeric(x), 3))
  
  return(all.stat)
}

# get t stats for sex
get_sex_stats <- function(model_list, factor_var = "sex") {
  
  stats_list <- lapply(model_list, function(mod) {
    
    # Estimated marginal means
    em <- emmeans(mod, specs = factor_var)
    
    em_df <- as.data.frame(em) %>%
      select(all_of(factor_var), emmean, SE, df) %>%
      rename(mean = emmean)
    
    # Pairwise contrast (F - M)
    contrast_obj <- contrast(em, method = "pairwise")
    contrast_df <- as.data.frame(contrast_obj) %>%
      select(contrast, estimate, SE, df, t.ratio, p.value)
    
    # Cohen's d for the contrast
    s <- sigma(mod)
    contrast_df$cohen_d <- contrast_df$estimate / s
    
    # Means in wide format
    em_wide <- em_df %>%
      pivot_wider(
        names_from  = all_of(factor_var),
        values_from = c(mean, SE),
        names_sep   = "_"
      )
    
    bind_cols(
      em_wide %>% select(-df),   # drop df from marginal means
      contrast_df
    )
  })
  
  # Add model names
  if (!is.null(names(model_list))) {
    stats_list <- Map(function(df, nm) {
      df$model <- nm
      df
    }, stats_list, names(model_list))
  }
  
  all_stats <- bind_rows(stats_list)
  
  # FDR correction on contrast p-values
  all_stats$p.fdr <- p.adjust(all_stats$p.value, method = "fdr")
  
  # Significance labels
  all_stats <- all_stats %>%
    mutate(
      sig = case_when(
        p.fdr < 0.05 ~ "YES",
        p.value < 0.05 ~ "yes",
        TRUE ~ NA_character_
      )
    )
  
  # Round numeric columns
  num_cols <- names(all_stats)[sapply(all_stats, is.numeric)]
  all_stats[num_cols] <- lapply(all_stats[num_cols], function(x) round(x, 3))
  
  all_stats
}




###############################################
#####         Filter significance        ######
###############################################

get.sig <- function(model.stat) {
  # filter rows where sig is not NA
  sig_rows <- model.stat[!is.na(model.stat$sig) & model.stat$sig %in% c("yes", "YES"), ]
  
  # return both region and term as a data.frame
  sig_info <- sig_rows[, c("model_name", "DV", "term", "p.value", "p.fdr", "sig")]
  
  return(sig_info)
}


###############################################
#####          Check Concurvity          ######
###############################################

check_and_filter_concurvity <- function(all_models, threshold = 0.8) {
  
  discarded_names <- list()
  
  for (list_name in names(all_models)) {
    
    cat("Checking", list_name, "\n")
    model_list <- all_models[[list_name]]
    
    high_count <- 0
    max_concur <- 0
    bad_model_names <- c()   # store only names
    
    for (model_name in names(model_list)) {
      mod <- model_list[[model_name]]
      
      cc <- concurvity(mod, full = TRUE)
      numeric_vals <- unlist(cc)
      numeric_vals <- numeric_vals[is.numeric(numeric_vals)]
      max_val <- if (length(numeric_vals) > 0) max(numeric_vals, na.rm = TRUE) else NA
      
      if (!is.na(max_val) && max_val > threshold) {
        bad_model_names <- c(bad_model_names, model_name)
        high_count <- high_count + 1
        if (max_val > max_concur) max_concur <- max_val
      }
    }
    
    # Summary print
    if (high_count > 0) {
      cat(sprintf(
        "In model list '%s': %d model(s) above threshold, max concurvity = %.3f\n",
        list_name, high_count, max_concur
      ))
      discarded_names[[list_name]] <- bad_model_names
    } else {
      cat(sprintf("No models above threshold in '%s'\n", list_name))
    }
  }
  
  cat("\nConcurvity check complete. Returning discarded model names.\n")
  return(discarded_names)
}


###############################################
#####           Get Derivatives          ######
###############################################


get.derivs <- function(modellist, labels, all_terms = FALSE) {
  result_list <- list()
  
  for (i in seq_along(labels)) {
    mod <- modellist[[i]]
    DV <- all.vars(formula(mod))[1]
    
    # Get all smooth terms
    all_smooth_terms <- rownames(summary(mod)$s.table)
    smooth_terms <- all_smooth_terms[startsWith(all_smooth_terms, "s")]
    
    # Determine n_terms dynamically and select the terms
    if (all_terms) {
      sm_terms <- smooth_terms
    } else {
      n_terms_use <- if (any(grepl(":sex", smooth_terms))) 2 else 1
      sm_terms <- head(smooth_terms, n_terms_use)
    }
    
    # Initialize deriv_vals with DV
    deriv_vals <- c(DV = DV)
    
    # Safely compute derivatives for all terms
    derivs_all <- tryCatch({
      lapply(sm_terms, function(sm_term) {
        derivs <- derivatives(mod, select = sm_term, type = "central")
        mean_deriv <- round(mean(derivs$.derivative, na.rm = TRUE), 3)
        clean_term <- gsub("[^A-Za-z]", "", sm_term)
        setNames(mean_deriv, paste0("deriv_", clean_term))
      })
    }, error = function(e) {
      message(sprintf("Skipping model '%s' due to error: %s — filling with NA", DV, e$message))
      # Return NA for all terms
      na_list <- setNames(rep(NA_real_, length(sm_terms)),
                          paste0("deriv_", gsub("[^A-Za-z]", "", sm_terms)))
      return(na_list)
    })
    
    # Combine the list of named derivatives into a single row
    deriv_vals <- c(DV = DV, unlist(derivs_all))
    
    # Convert to one-row data frame
    result_list[[i]] <- as.data.frame(as.list(deriv_vals), stringsAsFactors = FALSE)
  }
  
  # Combine all rows
  result_df <- bind_rows(result_list)
  
  # Convert numeric columns properly (everything except DV)
  result_df <- result_df %>%
    mutate(across(-DV, as.numeric))
  
  return(result_df)
}


