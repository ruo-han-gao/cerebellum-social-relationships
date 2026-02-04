library(rstudioapi)    # for getting script path
library(tidyverse)     # includes dplyr, ggplot2, tidyr, etc.
library(reshape2)      # for melt()
library(car)           # for vif()
library(interactions)  # for interact_plot()
library(pals)          # for coolwarm color

rm(list = ls())
options(scipen = 999, digits = 3)

codedir = dirname(getActiveDocumentContext()$path)
workdir = "/data/tu_gaor/scripts"
homedir = "/Users/ruohan/My_Repos/Cerebellum"
cwd_dir = homedir

# set the data set and DV of interest here
dataset_name = "hcpd"
yeo_of_interest <- "yeo7"

if (dataset_name == "hcpa") {
  sample_size <- 590
} else if (dataset_name == "hcpd") {
  sample_size <- 443
} else {
  print("dataset not recognized")
}

setwd(codedir)
data_dir = file.path(cwd_dir, "dat", dataset_name)
plot_dir = paste0(cwd_dir, "/plotting/assoc_plots/", yeo_of_interest)


# =============================================
#              Read in data                ####
# =============================================

dat <- read.csv(file.path(data_dir, paste0(dataset_name, "_brain_social_factors_n", sample_size, ".csv")))

# get all labels
groups <- list(
  fusion  = 5:36,
  yeo7 = 37:50,
  yeo17 = 51:84,
  social = 107:117
)

# get all labels
labels <- lapply(groups, function(cols) colnames(dat[, cols]))
fusion.labels = labels$fusion

yeo.labels <- labels[[yeo_of_interest]]

# rename yeo for better plotting
rename_yeo<- function(x) {
  hemi <- ifelse(grepl("^lh", x, ignore.case = TRUE), "L", "R")
  num  <- sub(".*_(\\d+)$", "\\1", x)
  paste0("N", num, hemi)
}
names(dat)[names(dat) %in% yeo.labels] <- rename_yeo(yeo.labels)

# update label names
labels <- lapply(groups, function(cols) colnames(dat[, cols]))
yeo.labels <- labels[[yeo_of_interest]]

pairs = expand.grid(
  cereb = fusion.labels,
  cortex = yeo.labels,
  stringsAsFactors = FALSE
)


# =============================================
#    Cerebellar cerebral correlations      ####
# =============================================

fusion_data = dat[,fusion.labels]
yeo_data = dat[,yeo.labels]

cor_matrix <- cor(fusion_data, yeo_data, 
                  use = "pairwise.complete.obs", 
                  method = "spearman")
range(cor_matrix) 
# 0.0262 0.5698 in hcpd
# 0.076 0.500 in hcpa

cor_long <- melt(cor_matrix)
cor_long <- cor_long %>%
  mutate(value_scaled = value / max(value, na.rm = TRUE))  # scale to 0-1

ggplot(cor_long) +
  geom_rect(aes(xmin = as.numeric(Var1) - 0.5*value_scaled,
                xmax = as.numeric(Var1) + 0.5*value_scaled,
                ymin = as.numeric(Var2) - 0.5*value_scaled,
                ymax = as.numeric(Var2) + 0.5*value_scaled,
                fill = value),
            color = NA) +
  scale_fill_gradientn(colors = pals::coolwarm(256),
                       limit = c(-0.6, 0.6), space = "Lab", breaks = c(-0.6, 0, 0.6),
                       name = expression(italic(r))) +
  scale_x_continuous(breaks = 1:length(unique(cor_long$Var1)), 
                     labels = unique(cor_long$Var1)) +
  scale_y_continuous(breaks = 1:length(unique(cor_long$Var2)), 
                     labels = unique(cor_long$Var2)) +
  coord_fixed() +
  theme_minimal() +
  theme(
    text = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "", y = "", title = "Raw Correlation Heatmap")

ggsave(file.path(plot_dir, paste0(dataset_name, "_cc_raw_correlation_", yeo_of_interest, ".png")), width=8, height=6, dpi=300)

# =============================================
#              Check normality             ####
# =============================================

# function
check_normality <- function(data, margin = NULL) {
  # margin = 2 → check columns
  # margin = 1 → check rows
  
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data.frame or matrix.")
  }
  
  # convert to matrix for consistency
  mat <- as.matrix(data)
  
  results <- lapply(1:dim(mat)[margin], function(i) {
    x <- if (margin == 1) as.numeric(mat[i, ]) else as.numeric(mat[, i])
    test <- shapiro.test(x)
    data.frame(
      Target = if (margin == 1) rownames(mat)[i] else colnames(mat)[i],
      W = test$statistic,
      p_value = test$p.value,
      Normal = ifelse(test$p.value < 0.05, "no", "yes")
    )
  })
  
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  colnames(results_df)[1] = "Parcel"
  
  return(results_df)
}

# check normality
normal.fusion <- check_normality(fusion_data, margin = 2)
table(normal.fusion$Normal)
normal.yeo <- check_normality(yeo_data, margin = 2)
table(normal.yeo$Normal)

# normal distribution not met in majority. However, giving the large sample size, I will ignore it and continue.


# =============================================
#         Regression Preparation           ####
# =============================================

# contrast for categorical variables
dat$sex <- as.factor(dat$sex)
contrasts(dat$sex) <- contr.sum(2)
colnames(contrasts(dat$sex)) = c("F_M")
contrasts(dat$sex) # grand mean across F and M as reference
# positive coef: females > males; negative beta: males > females


# to scale continuous variables for standardized beta
contin_vars <- c(
  "icv", "age", "sr", fusion.labels, yeo.labels
)

# Also centering continuous predictors to reduce multicollinearity in interactions 
dat[, contin_vars] <- scale(dat[, contin_vars], center = TRUE, scale = TRUE)

# =============================================
#             Run regression models        ####
# =============================================

results_list <- lapply(seq_len(nrow(pairs)), function(i) {
  
  var1 <- pairs$cortex[i]   # DV
  var2 <- pairs$cereb[i]    # IV
  
  # Build formula
  fmla <- as.formula(paste0(var1, " ~ ", var2, " * sr + age + sex + icv"))
  
  # Fit model
  model <- lm(fmla, data = dat)
  
  # Extract rows of interest from summary
  coef_df <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    # filter(term %in% c("sr", var2, paste0(var2, ":sr"), "age", "sex1")) %>%
    mutate(DV = var1,
           IV = var2)
  
  # Return both the model and coefficient table
  list(DV = var1, model = model, coef = coef_df)
})

# Combine all coefficient tables into a single data frame
result_df <- bind_rows(lapply(results_list, function(x) x$coef))

# =============================================
#           Check multicollinearity        ####
# =============================================

# Preallocate a list to store VIF results
vif_results <- vector("list", length(results_list))

for (i in seq_along(results_list)) {
  model <- results_list[[i]]$model
  if (!is.null(model)) {
    # Include interactions but suppress repeated warning
    vif_results[[i]] <- suppressMessages(vif(model, type = "terms"))
  } else {
    vif_results[[i]] <- NA
  }
}

# Check if any VIF > 5
vif_over_5 <- sapply(vif_results, function(v) {
  if (is.numeric(v)) any(v > 5) else FALSE
})

if (any(vif_over_5)) {
  warning_models <- which(vif_over_5)
  message(
    "\n High VIF detected (VIF > 5) in models:\n",
    paste(warning_models, collapse = ", "),
    "\nCheck predictors for multicollinearity."
  )
} else {
  message("No VIF values exceed 5. Multicollinearity OK.")
}


# =============================================
#          Find significant results        ####
# =============================================

# apply fdr correction
result_stats <- result_df %>%
  mutate(
    p_fdr = p.adjust(`Pr(>|t|)`, method = "fdr"),
    sig = case_when(
      grepl(":sr", term) & p_fdr < 0.05        ~ "YES:int",   # interaction, FDR significant
      grepl(":sr", term) & `Pr(>|t|)` < 0.05   ~ "yes:int",   # interaction, uncorrected sig
      !grepl(":sr", term) & p_fdr < 0.05       ~ "YES",       # main effect, FDR sig
      !grepl(":sr", term) & `Pr(>|t|)` < 0.05  ~ "yes",       # main effect, uncorrected sig
      TRUE                                       ~ NA_character_
    )
  )

results_sig <- result_stats %>%
  filter(!is.na(sig) & tolower(sig) == "yes:int") %>%
  select(DV, term, Estimate, `t value`, `Pr(>|t|)`, p_fdr, sig)


write.csv(result_stats, row.names = FALSE, 
          file = file.path(data_dir, "cc_assoc", 
                           paste0(dataset_name, "_covar_fusion_", yeo_of_interest, ".csv")))

write.csv(results_sig, row.names = FALSE, 
          file = file.path(data_dir, "cc_assoc", 
                           paste0(dataset_name, "_sig_covar_fusion_", yeo_of_interest, ".csv")))

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

# heatmap 1
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
  scale_fill_gradientn(
    colors = pals::coolwarm(256), limits = c(-0.3, 0.3), breaks = c(-0.3, 0, 0.3),
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

ggsave(file.path(plot_dir, paste0(dataset_name, "_heatmap_fusion_sr_int_", yeo_of_interest, ".png")), width=8, height=6, dpi=300)

# heatmap 1
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

