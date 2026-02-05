library(rstudioapi)    # for getting script path
library(tidyverse)     # includes dplyr, ggplot2, tidyr, etc.
library(reshape2)      # for melt()
library(car)           # for vif()
library(interactions)  # for interact_plot()
library(pals)          # for coolwarm color

source("1_covar_analysis.R")

# =============================================
#        Plot main effect term -SR        ####
# =============================================

fusion_terms <- result_stats %>%
  filter(term == "sr")

# Prepare data
heatmap_data <- fusion_terms %>%
  mutate(beta = ifelse(`Pr(>|t|)` < 0.05, `Estimate`, NA)) %>%
  select(DV, IV, beta, p_fdr) %>%  # keep p_fdr for stars
  mutate(
    DV = factor(DV, levels = unique(DV)),
    IV = factor(IV, levels = unique(IV))
  )
range(heatmap_data$beta, na.rm = TRUE)

# heatmap
ggplot(heatmap_data, aes(x = IV, y = DV, fill = beta)) +
  # Base layer: all tiles with light gray border
  geom_tile(color = "lightgray") +
  # Overlay: add a white "*" for FDR-significant tiles
  geom_text(
    data = filter(heatmap_data, p_fdr < 0.05),
    aes(label = "*"),
    color = "white",
    size = 5, hjust = 0.5, vjust = 0.8
  ) +
  scale_fill_gradientn(colors = pals::coolwarm(256),
    limits = c(-0.1, 0.1), breaks = c(-0.1, 0, 0.1),
    na.value = "white"
  ) +
  labs(x = "", y = "", title = toupper(paste0(dataset_name)), fill = "\u03B2"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_text(size = 12, family = "Arial", face = "italic"),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(plot_dir, paste0(dataset_name, "_cc_assoc_fusion_sr_on", yeo_of_interest, ".png")), width=8, height=6, dpi=300)

# =============================================
#            Plot main effect Cerebel      ####
# =============================================

fusion_terms <- result_stats %>%
  filter(grepl("[1-9]", term) & !grepl(":", term))

# Prepare data
heatmap_data <- fusion_terms %>%
  mutate(beta = ifelse(`Pr(>|t|)` < 0.05, `Estimate`, NA)) %>%
  select(DV, IV, beta, p_fdr) %>%  # keep p_fdr for stars
  mutate(
    DV = factor(DV, levels = unique(DV)),
    IV = factor(IV, levels = unique(IV))
  )
range(heatmap_data$beta, na.rm = TRUE)

# heatmap
ggplot(heatmap_data, aes(x = IV, y = DV, fill = beta)) +
  # Base layer: all tiles with light gray border
  geom_tile(color = "lightgray") +
  # Overlay: add a white "*" for FDR-significant tiles
  geom_text(
    data = filter(heatmap_data, p_fdr < 0.05),
    aes(label = "*"),
    color = "white",
    size = 5, hjust = 0.5, vjust = 0.8
  ) +
  scale_fill_gradientn(colors = pals::coolwarm(256),
                       limits = c(-0.3, 0.3), breaks = c(-0.3, 0, 0.3),
                       na.value = "white"
  ) +
  labs(x = "", y = "", title = toupper(paste0(dataset_name)), fill = "\u03B2"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_text(size = 12, family = "Arial", face = "italic"),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(plot_dir, paste0(dataset_name, "_cc_assoc_fusion_on", yeo_of_interest, ".png")), width=8, height=6, dpi=300)

# =============================================
#           Plot main effect - sex         ####
# =============================================

sex_term <- result_stats %>%
  filter(term == "sexF_M")

# Prepare data
heatmap_data <- sex_term %>%
  mutate(beta = ifelse(`Pr(>|t|)` < 0.05, `Estimate`, NA)) %>%
  select(DV, IV, beta, p_fdr) %>%  # keep p_fdr for stars
  mutate(
    DV = factor(DV, levels = unique(DV)),
    IV = factor(IV, levels = unique(IV))
  )
range(heatmap_data$beta, na.rm = TRUE)

# heatmap
ggplot(heatmap_data, aes(x = IV, y = DV, fill = beta)) +
  # Base layer: all tiles with light gray border
  geom_tile(color = "lightgray") +
  # Overlay: add a white "*" for FDR-significant tiles
  geom_text(
    data = filter(heatmap_data, p_fdr < 0.05),
    aes(label = "*"),
    color = "white",
    size = 5, hjust = 0.5, vjust = 0.8
  ) +
  scale_fill_gradientn(colors = pals::coolwarm(256),
                       limits = c(-0.3, 0.3), breaks = c(-0.3, 0, 0.3),
                       na.value = "white"
  ) +
  labs(x = "", y = "", title = toupper(paste0(dataset_name)), fill = "\u03B2"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_text(size = 12, family = "Arial", face = "italic"),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(plot_dir, paste0(dataset_name, "_cc_assoc_sex_", yeo_of_interest, ".png")), width=8, height=6, dpi=300)


# =============================================
#           Plot main effect - age         ####
# =============================================

age_term <- result_stats %>%
  filter(term == "age")

# Prepare data
heatmap_data <- age_term %>%
  mutate(beta = ifelse(`Pr(>|t|)` < 0.05, `Estimate`, NA)) %>%
  select(DV, IV, beta, p_fdr) %>%  # keep p_fdr for stars
  mutate(
    DV = factor(DV, levels = unique(DV)),
    IV = factor(IV, levels = unique(IV))
  )
range(heatmap_data$beta, na.rm = TRUE)

# heatmap
ggplot(heatmap_data, aes(x = IV, y = DV, fill = beta)) +
  # Base layer: all tiles with light gray border
  geom_tile(color = "lightgray") +
  # Overlay: add a white "*" for FDR-significant tiles
  geom_text(
    data = filter(heatmap_data, p_fdr < 0.05),
    aes(label = "*"),
    color = "white",
    size = 5, hjust = 0.5, vjust = 0.8
  ) +
  scale_fill_gradientn(colors = pals::coolwarm(256),
                       limits = c(-0.7, 0.7), breaks = c(-0.5, 0, 0.5),
                       na.value = "white"
  ) +
  labs(x = "", y = "", title = toupper(paste0(dataset_name)), fill = "\u03B2"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_text(size = 12, family = "Arial", face = "italic"),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(plot_dir, paste0(dataset_name, "_cc_assoc_age_", yeo_of_interest, ".png")), width=8, height=6, dpi=300)




# =============================================
#            Plot interaction term         ####
# =============================================

interaction_terms <- result_stats %>%
  filter(grepl(":sr$", term))

# Prepare data
heatmap_data <- interaction_terms %>%
  mutate(beta = ifelse(`Pr(>|t|)` < 0.05, `Estimate`, NA)) %>%
  select(DV, IV, beta, p_fdr) %>%  # keep p_fdr for black border
  mutate(
    DV = factor(DV, levels = unique(DV)),
    IV = factor(IV, levels = unique(IV))
  )
range(heatmap_data$beta, na.rm = TRUE)

# heatmap
ggplot(heatmap_data, aes(x = IV, y = DV, fill = beta)) +
  # Base layer: all tiles with light gray border
  geom_tile(color = "lightgray") +
  # Overlay: add a white "*" for FDR-significant tiles
  geom_text(
    data = filter(heatmap_data, p_fdr < 0.05),
    aes(label = "*"),
    color = "white",
    size = 5, hjust = 0.5, vjust = 0.8    # horizontal center  # adjust size as needed
  ) +
  scale_fill_gradientn(colors = pals::coolwarm(256),limits = c(-0.1, 0.1), breaks = c(-0.1, 0, 0.1),
    na.value = "white"
  ) +
  labs(x = "", y = "", title = toupper(paste0(dataset_name)), fill = "\u03B2") +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_text(size = 12, family = "Arial", face = "italic"),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(plot_dir, paste0(dataset_name, "_heatmap_fusion_sr_int_", yeo_of_interest, "_1.png")), width=8, height=6, dpi=300)



# =============================================
# Plot individual interaction effect       ####
# =============================================

individual_plot_data <- heatmap_data %>%
  filter(!is.na(beta))

# Loop over rows of individual_plot_data:
for (i in seq_len(nrow(individual_plot_data))) {
  
  dv <- as.character(individual_plot_data$DV[i])     # e.g., "rh_17Networks_15"
  iv  <- as.character(individual_plot_data$IV[i])     # e.g., "D3L"
  fdr_flag <- isTRUE(individual_plot_data$p_fdr[i] < 0.05)
  
  index = which(pairs$cereb == iv & pairs$cortex == dv)
  # Get the model for this pair
  model <- results_list[[index]]$model
  
  # Build interaction plot
  p <- interact_plot(
    model,
    pred = !!sym(iv),      
    modx = sr,              # moderator = sr
    plot.points = TRUE,
    partial.residuals = TRUE
  )
  
  for (a in c("colour", "linetype")) {
    sc <- p$scales$get_scales(a)
    if (!is.null(sc) && !is.null(sc$name)) {
      sc$name <- toupper(sc$name)
    }
  }
  
  p <- p +
    labs(x = iv, y = dv, title = toupper(dataset_name)) +
    theme(
      text = element_text(size = 18),
      axis.text = element_text(size = 16)
    )
  
  if (fdr_flag) {
    p <- p +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = "*",
        hjust = 1.2, vjust = 1.2,
        size = 11,
        color = "black"
      )
  }
  
  # Output filename
  out_file <- file.path(
    plot_dir,
    paste0(dataset_name, "_", iv, "_",yeo_of_interest, "_", dv, "_interaction.png")
  )
  
  print(paste(
    "Saving plot", i, "/", nrow(individual_plot_data),
    ": IV =", iv, ", DV =", dv, "for", dataset_name, yeo_of_interest))
  
  
  # Save plot without rescaling content
  ggsave(
    filename = out_file,
    plot = p,
    width = 5,
    height = 4,
    dpi = 300,
    scaling = 1
  )
}

