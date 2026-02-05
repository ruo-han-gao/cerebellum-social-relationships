library(dplyr)
library(purrr)
library(RColorBrewer)
library(colorspace)
library(rstudioapi) # gets path of script

codedir = dirname(getActiveDocumentContext()$path)
setwd(codedir)

# if the gam model results are not in the environment, then to plot, first run
source("2_gam_for_social.R")


plot_dir = file.path(cwd_dir, dataset_name, "from_gam/plots/curves")

# =============================================
# Set plotting colors                      ####
# =============================================

parcel_colors_clean <- c(
  A1 = "#898AFD",
  A2 = "#544EE3",
  A3 = "#60FFFE",
  
  M1 = "#804EA1",
  M2 = "#9D29E0",
  M3 = "#65B5F8",
  M4 = "#62FF74",
  
  D1 = "#AB3F48",
  D2 = "#FE6BB5",
  D3 = "#F054FF",
  D4 = "#D01E82",
  
  S1 = "#FDC122",
  S2 = "#DBFF48",
  S3 = "#E3DE64",
  S4 = "#E4D8AC",
  S5 = "#53BD60"
)

parcel_colors <- c(
  setNames(parcel_colors_clean, paste0(names(parcel_colors_clean), "L")),
  setNames(parcel_colors_clean, paste0(names(parcel_colors_clean), "R"))
)
# =============================================
# Function to get plotting data            ####
# =============================================

get_gam_plot_data <- function(model_list, label_of_interest, term_to_plot) {
  
  names(model_list) <- label_of_interest
  
  # Get stats (keep ALL rows)
  stats <- get.gam.stats(model_list)
  
  # Lookup table for significance for this term
  sig_lookup <- stats %>%
    filter(term == term_to_plot) %>%
    select(DV, sig) %>%
    distinct()
  
  get_sig <- function(dv) {
    s <- sig_lookup$sig[sig_lookup$DV == dv]
    if (length(s) == 0) NA else s[[1]]
  }
  
  smooth_list <- list()
  residual_list <- list()
  
  for (mod_name in names(model_list)) {
    model <- model_list[[mod_name]]
    
    sig_val <- get_sig(mod_name)
    is_sig  <- !is.na(sig_val)
    is_fdr  <- sig_val %in% c("Y", "YES")
    
    # Smooth curve data
    pdat <- plot(model, pages = 1, seWithMean = TRUE)
    
    smooth_df <- map_dfr(seq_along(pdat), function(i) {
      data.frame(
        term   = gsub(",[0-9.]+", "", pdat[[i]]$ylab),
        x      = pdat[[i]]$x,
        fit    = pdat[[i]]$fit / 1000,
        se     = pdat[[i]]$se  / 1000,
        Parcel = mod_name,
        sig    = sig_val,
        is_sig = is_sig,
        is_fdr = is_fdr
      )
    }) %>%
      filter(term == term_to_plot)
    
    smooth_list[[mod_name]] <- smooth_df
    
    # Partial residuals
    term_names <- colnames(predict(model, type = "terms"))
    if (term_to_plot %in% term_names) {
      
      xvar <- gsub("s\\((.*)\\).*", "\\1", term_to_plot)
      
      resid_df <- data.frame(
        x = model$model[[xvar]],
        resid_partial = (residuals(model, type = "response") +
                           predict(model, type = "terms")[, term_to_plot]) / 1000,
        Parcel = mod_name,
        sig    = sig_val,
        is_sig = is_sig,
        is_fdr = is_fdr
      )
      
      residual_list[[mod_name]] <- resid_df
    }
  }
  
  all_smooth    <- bind_rows(smooth_list)
  all_residuals <- bind_rows(residual_list)
  
  list(smooth = all_smooth, residuals = all_residuals)
}


# =============================================
# Plot curves                              ####
# =============================================

residual <- TRUE


# =============================================
# Plot selected curves for age            ####
# =============================================

model_list = final_models$sr.mod.srsmoo
terms = c("s(age)")

hcpd_age_select = c("A2R","M1L", "D3R", "S2L")
hcpa_age_select = c("A2R","M3R", "D1L", "S5R")

for (term in terms) {
  plot_data <- get_gam_plot_data(
    model_list = model_list,
    label_of_interest = label_of_interest,
    term_to_plot = term
  )
  selected_parcels_name <- paste0(dataset_name, "_age_select")
  selected_parcels <- get(selected_parcels_name)
  
  smooth_data <- plot_data$smooth
  resid_data <- plot_data$residuals
  
  smooth_data_sub <- plot_data$smooth %>%
    filter(Parcel %in% selected_parcels)
  
  resid_data_sub <- plot_data$residuals %>%
    filter(Parcel %in% selected_parcels)
  
  # Create the plot and assign to a variable
  p <- ggplot() +
    geom_point(data = resid_data_sub, aes(x = x, y = resid_partial, color = Parcel), alpha = 0.2, size = 0.1) +
    geom_line(data = smooth_data_sub, aes(x = x, y = fit, color = Parcel), linewidth = 0.7) +
    scale_color_manual(values = parcel_colors) +
    theme_minimal() +
    labs(
      x = "Age (centered)",
      y = expression(Volume~(cm^3)),
      title = paste("")
    ) +
    theme(
      axis.title = element_text(size = 11, color = "black", family = "Arial"),
      axis.text  = element_text(size = 11, color = "black", family = "Arial"),
      legend.text  = element_text(size = 11, family = "Arial"),
      legend.title = element_text(size = 11, family = "Arial"),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  
  # Save the plot to file
  ggsave(
    filename = file.path(plot_dir, paste0("GAM_curve_", dataset_name, "_", gsub("[^A-Za-z0-9]", "_", term), ".png")),
    plot = p, width = 3.5, height = 2.5, dpi = 300)
}

# =============================================
# Individual curves for age significant    ####
# =============================================

# no residuals 
for (term in terms) {
  
  plot_data <- get_gam_plot_data(
    model_list = model_list,
    label_of_interest = label_of_interest,
    term_to_plot = term
  )
  
  smooth_data <- plot_data$smooth
  
  parcels <- sort(unique(smooth_data$Parcel))
  
  for (parcel in parcels) {
    
    smooth_one <- smooth_data %>% filter(Parcel == parcel)
    
    if (nrow(smooth_one) == 0) next
    
    # Significant? (treat NA as FALSE)
    sig_flag <- any(smooth_one$is_sig %in% TRUE)
    fdr_flag <- any(smooth_one$is_fdr %in% TRUE)
    
    lt <- if (isTRUE(sig_flag)) "solid" else "dashed"
    
    col <- parcel_colors[[parcel]]
    if (is.null(col)) col <- "black"
    
    # star position (top-right, inside panel)
    star_df <- if (isTRUE(fdr_flag)) {
      data.frame(
        x = max(smooth_one$x, na.rm = TRUE),
        y = max(smooth_one$fit + 2 * smooth_one$se, na.rm = TRUE),
        label = "*"
      )
    } else NULL
    
    p <- ggplot() +
      geom_ribbon(
        data = smooth_one,
        aes(x = x, ymin = fit - 2 * se, ymax = fit + 2 * se),
        fill = col,
        alpha = 0.18
      ) +
      geom_line(
        data = smooth_one,
        aes(x = x, y = fit),
        linewidth = 0.7,
        color = col,
        linetype = lt
      ) +
      { if (!is.null(star_df))
        geom_text(
          data = star_df,
          aes(x = x, y = y, label = label),
          size = 6,
          vjust = 1,
          hjust = 1,
          color = "black"
        )
      } +
      theme_minimal() +
      labs(
        x = "Age (centered)",
        y = expression(Volume~(cm^3)),
        title = parcel
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 11, color = "black", family = "Arial"),
        axis.text  = element_text(size = 11, color = "black", family = "Arial"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)
      )
    
    ggsave(
      filename = file.path(
        plot_dir,
        paste0(
          "GAM_curve_", dataset_name, "_",
          gsub("[^A-Za-z0-9]", "_", term), "_",
          gsub("[^A-Za-z0-9]", "_", parcel), ".png"
        )
      ),
      plot = p, width = 3.5, height = 2.5, dpi = 300
    )
  }
}


# =============================================
# Plot curve for social scores             ####
# =============================================

get_xlab <- function(term) {
  if (grepl("sr",  term, ignore.case = TRUE))  return("SR")
  if (grepl("rej", term, ignore.case = TRUE))  return("Rejection")
  if (grepl("hos", term, ignore.case = TRUE))  return("Hostility")
  if (grepl("fri", term, ignore.case = TRUE))  return("Friendship")
  if (grepl("lon", term, ignore.case = TRUE))  return("Loneliness")
  return("")  # fallback
}


plot_gam_terms_sig <- function(terms, model_list) {
  
  for (term in terms) {
    
    plot_data <- get_gam_plot_data(
      model_list = model_list,
      label_of_interest = label_of_interest,
      term_to_plot = term
    )
    
    smooth_data <- plot_data$smooth %>% filter(is_sig %in% TRUE)
    resid_data  <- plot_data$residuals %>% filter(is_sig %in% TRUE)
    
    p <- ggplot() +
      geom_point(
        data = resid_data,
        aes(x = x, y = resid_partial, color = Parcel),
        alpha = 0.1, size = 0.1
      ) +
      geom_line(
        data = smooth_data,
        aes(x = x, y = fit, color = Parcel),
        linewidth = 0.7
      ) +
      scale_color_manual(values = parcel_colors) +
      theme_minimal() +
      labs(
        x = get_xlab(term),
        y = expression(Volume~(cm^3)),
        title = ""
      )
    
    ggsave(
      filename = file.path(
        plot_dir,
        paste0(
          "GAM_curve_", dataset_name, "_",
          gsub("[^A-Za-z0-9]", "_", term), ".png"
        )
      ),
      plot = p,
      width = 3,
      height = 2,
      dpi = 300
    )
  }
}

# =============================================
# Plot curves for social relationship.     ####
# =============================================
model_list = final_models$sr.mod.srsmoo
terms = c("s(sr)")

# plot general SR
plot_gam_terms_sig(terms, model_list)


# Plot curves for social score:sex
model_list <- final_models$sr.mod.srsex
terms <- c("s(sr):sexM", "s(sr):sexF")

plot_gam_terms_sig(terms, model_list)


# =============================================
# Plot curves for raw_loneli               ####
# =============================================

model_list = final_models$lon.mod.lonsmoo
terms = c("s(raw_loneli)")

# plot general raw_loneli effect
plot_gam_terms_sig(terms, model_list)


# Plot curves for raw_loneli:sex
model_list <- final_models$lon.mod.lonsex
terms <- c("s(raw_loneli):sexF", "s(raw_loneli):sexM")

plot_gam_terms_sig(terms, model_list)

# =============================================
# Plot curves for rejection                ####
# =============================================

model_list = final_models$rej.mod.rejsmoo
terms = c("s(raw_rej)")

# plot general rejection effect
plot_gam_terms_sig(terms, model_list)


# Plot curves for rejection:sex
model_list <- final_models$rej.mod.rejsex
terms <- c("s(raw_rej):sexF", "s(raw_rej):sexM")

plot_gam_terms_sig(terms, model_list)

# =============================================
# Plot curves for friendship               ####
# =============================================

model_list = final_models$fri.mod.frismoo
terms = c("s(raw_friends)")

# plot general friendship effect
plot_gam_terms_sig(terms, model_list)


# Plot curves for friendship:sex
model_list <- final_models$fri.mod.frisex
terms <- c("s(raw_friends):sexF", "s(raw_friends):sexM")

plot_gam_terms_sig(terms, model_list)

# =============================================
# Plot curves for raw_hostility            ####
# =============================================

model_list = final_models$hos.mod.hossmoo
terms = c("s(raw_hostility)")

# plot general friendship effect
plot_gam_terms_sig(terms, model_list)


# Plot curves for friendship:sex
model_list <- final_models$hos.mod.hossex
terms <- c("s(raw_hostility):sexF", "s(raw_hostility):sexM")

plot_gam_terms_sig(terms, model_list)


