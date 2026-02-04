library(rstudioapi) # gets path of script
library(mgcv)
library(dplyr)
library(purrr) # for combine multiple dfs (reduce)
library(tidyr)
library(broom)
library(ggplot2)
library(gratia) # derivatives
library(emmeans) # for estimated means and stats


# clear the environment
rm(list = ls())
options(scipen = 999, digits = 3)

codedir = dirname(getActiveDocumentContext()$path)
workdir = "/data/tu_gaor/scripts/dat"
homedir = "/Users/ruohan/My_Repos/Cerebellum/dat"
cwd_dir = homedir

# set the data set and DV of interest here
dataset_name = "hcpd"
DV_of_interest <- "fusion"

if (dataset_name == "hcpa") {
  sample_size <- 590
} else if (dataset_name == "hcpd") {
  sample_size <- 443
} else {
  print("dataset not recognized")
}

data_dir = file.path(cwd_dir, dataset_name)
setwd(codedir)


# To run this script, first run `gam_functions.R` for necessary functions
source("2_gam_functions.R")

# =============================================
#            Read in data                  ####
# =============================================

# get cerebellar and cortical data
dat <- read.csv(file.path(data_dir, paste0(dataset_name, "_brain_social_factors_n", sample_size, ".csv")))
mean(dat$age) # M.hcpd  = 14.9; M.hcpa = 58
sd(dat$age)   # SD.hcpd = 3.93; SD.hcpa = 14.1

# get all labels
groups <- list(
  fusion  = 5:36,
  yeo7 = 37:50,
  yeo17 = 51:84,
  social = 107:117
)

labels <- lapply(groups, function(cols) colnames(dat[, cols]))
social.labels = labels$social

label_of_interest <- labels[[DV_of_interest]]

# =============================================
#         Regression Preparation           ####
# =============================================

# contrast for categorical variables
dat$sex <- as.factor(dat$sex)
contrasts(dat$sex) <- contr.sum(2)
colnames(contrasts(dat$sex)) = c("F_M")
contrasts(dat$sex) # grand mean across F and M as reference
# positive coef: females > males; negative beta: males > females

# Define continuous covariates
concovs <- c(
  "icv", "age", "raw_friends", "raw_emosup", "raw_loneli", "raw_rej",
  "raw_hostility", "sr", "Loneliness", "Hostility", "EmoSupt",
  "Friends", "Rejection"
)
mean(dat$age)
# Centering continuous predictors to reduce multicollinearity in interactions 
dat[, concovs] <- scale(dat[, concovs], center = TRUE, scale = FALSE)

# Check correlation between continuous predictors
cor_age_icv <- cor.test(dat$icv, dat$age)
print(cor_age_icv)
# hcpd: p-value = 0.6339, cor =0.01149961 → no multicollinearity concern
# hcpa: p-value = 0.0000, cor =0.213 →  multicollinearity concern, check concurvity 

cor_age_social <- data.frame(
  r = sapply(social.labels, function(f)
    cor(dat[[f]], dat[["age"]], use = "complete.obs")),
  p = sapply(social.labels, function(f)
    cor.test(dat[[f]], dat[["age"]])$p.value)
)
print(cor_age_social)

# =============================================
# Can age and sex predict social score?    ####
# =============================================

mod.social = fit.models(formula = "s(age) + sex", labels = social.labels)
stat_social = get.gam.stats(mod.social)
plot(mod.social[[1]], pages = 1, residuals = TRUE, pch = 1, cex = 0.5) # raw_friends
plot(mod.social[[2]], pages = 1, residuals = TRUE, pch = 1, cex = 0.5) # raw_emosup
plot(mod.social[[3]], pages = 1, residuals = TRUE, pch = 1, cex = 0.5) # raw_loneli
plot(mod.social[[4]], pages = 1, residuals = TRUE, pch = 1, cex = 0.5) # raw_rej
plot(mod.social[[5]], pages = 1, residuals = TRUE, pch = 1, cex = 0.5) # raw_hostility
plot(mod.social[[6]], pages = 1, residuals = TRUE, pch = 1, cex = 0.5) # sr

# =============================================
# Fit and compare covarience models        ####
# =============================================

## icv ====
mod.icv = fit.models(formula = "icv", labels = label_of_interest)
mod.icv.smoo = fit.models(formula = "s(icv)", labels = label_of_interest)
compare.icvsmoo = compare.models(mod.icv, mod.icv.smoo)
icv.smoo.better = better.models(compare.icvsmoo)

## sex ====
mod.sex = fit.models(formula = "s(icv) + sex", labels = label_of_interest)
compare.sex = compare.models(mod.icv.smoo, mod.sex)
sex.better = better.models(compare.sex)

## age ====
mod.age = fit.models(formula = "age + s(icv) + sex", labels = label_of_interest)
compare.age = compare.models(mod.sex, mod.age)
age.better = better.models(compare.age)

## age smooth ====
mod.agesmoo = fit.models(formula = "s(age) + s(icv) + sex", labels = label_of_interest)
compare.agesmoo = compare.models(mod.sex, mod.agesmoo)
agesmoo.better = better.models(compare.agesmoo)

## age sex interaction ====
mod.agesex = fit.models(formula = "s(age, by = sex) + s(icv) + sex", labels = label_of_interest)
compare.agesex = compare.models(mod.sex, mod.agesmoo)
agesex.better = better.models(compare.agesex)


# =============================================
# Fit and compare social var models        ####
# =============================================

sr_models = fit_hirac_models("sr", mod.agesmoo, label_of_interest)

range(sr_models$comparisons$compare.srsmoo$delta_R2)
# hcpd: -0.010  0.013;   hcpa: -0.001  0.003
range(sr_models$comparisons$compare.srsex$delta_R2)
# hcpd: -0.003  0.021;   hcpa: -0.003  0.008
range(sr_models$comparisons$compare.srage1$delta_R2)
# hcpd: -0.012  0.013;   hcpa: -0.003  0.003
range(sr_models$comparisons$compare.srage2$delta_R2)
# hcpd: 0.000 0.017;     hcpa: 0.000 0.012

lon_models = fit_hirac_models("raw_loneli", mod.agesmoo, label_of_interest)
rej_models = fit_hirac_models("raw_rej", mod.agesmoo, label_of_interest)
hos_models = fit_hirac_models("raw_hostility", mod.agesmoo, label_of_interest)
emo_models = fit_hirac_models("raw_emosup", mod.agesmoo, label_of_interest)
fri_models = fit_hirac_models("raw_friends", mod.agesmoo, label_of_interest)

# lon_models = fit_hirac_models("Loneliness", mod.agesmoo, label_of_interest)
# rej_models = fit_hirac_models("Rejection", mod.agesmoo, label_of_interest)
# hos_models = fit_hirac_models("Hostility", mod.agesmoo, label_of_interest)
# emo_models = fit_hirac_models("EmoSupt", mod.agesmoo, label_of_interest)
# fri_models = fit_hirac_models("Friends", mod.agesmoo, label_of_interest)

# =============================================
# Filter out models with high concurvity   ####
# =============================================

all_models <- c(
  cov = list(mod.agesmoo),
  sr  = sr_models$models,
  lon = lon_models$models,
  rej = rej_models$models,
  hos = hos_models$models,
  emo = emo_models$models,
  fri = fri_models$models
)

# Check concurvity for all predictors
discarded <- check_and_filter_concurvity(all_models, threshold = 0.8)
discarded_names <- sub(".*\\.", "", names(discarded))
discard_pattern <- paste(discarded_names, collapse = "$|") |> paste0("$")

# =============================================
# Stats from final models                  ####
# =============================================

final_models <- all_models[!grepl(discard_pattern, names(all_models))]

# all stats
gam_stat_list <- lapply(names(final_models), function(name) {
  cat("Get GAM stats for", name, "\n")
  get.gam.stats(final_models[[name]], model_name = name)
})
gam_stats_all <- bind_rows(gam_stat_list)
gam_sig = get.sig(gam_stats_all)


# stats for linear interaction term
int_stat_list <- lapply(names(final_models), function(name) {
  cat("Get interation stats for", name, "\n")
  get.int.stats(final_models[[name]], model_name = name)
})
int_stat_all <- bind_rows(int_stat_list)

# t stats for sex
covar_model = final_models$cov
sex_stats = get_sex_stats(covar_model, factor_var = "sex")


write.csv(gam_stats_all, row.names = FALSE, 
          file = file.path(data_dir, "from_gam", 
                           paste0(dataset_name, "_GAM_stats_", DV_of_interest, ".csv")))

write.csv(gam_sig, row.names = FALSE, 
          file = file.path(data_dir, "from_gam", 
                           paste0(dataset_name, "_GAM_signif_", DV_of_interest, ".csv")))

write.csv(int_stat_all, row.names = FALSE, 
          file = file.path(data_dir, "from_gam", 
                           paste0(dataset_name, "_GAM_linear_int_term", DV_of_interest, ".csv")))

write.csv(sex_stats, row.names = FALSE, 
          file = file.path(data_dir, "from_gam", 
                           paste0(dataset_name, "_GAM_sex_stats_", DV_of_interest, ".csv")))

# =============================================
# Partial R2 from a reduced or added mode  ####
# =============================================

## without age ====
mod._age = fit.models(formula = "s(icv) + sex", labels = label_of_interest)
compare._age = compare.models(mod._age, mod.agesmoo)
worse._age = better.models(compare._age)

## without sex ====
mod._sex = fit.models(formula = "s(icv) + s(age)", labels = label_of_interest)
compare._sex = compare.models(mod._sex, mod.agesmoo)
worse._sex = better.models(compare._sex)

all_comparison = c(
  age = list(compare._age),
  sex = list(compare._sex),
  sr  = sr_models$comparisons,
  lon = lon_models$comparisons,
  rej = rej_models$comparisons,
  hos = hos_models$comparisons,
  emo = emo_models$comparisons,
  fri = fri_models$comparisons
)

final_compare <- all_comparison[!grepl(discard_pattern, names(all_comparison))]

# Initialize list to store R2 values
r2_list <- list()

# Keep the DV column from age as 'region'
region <- final_compare$age$DV

for (nm in names(final_compare)) {
  df <- final_compare[[nm]]
  
  r2_val <- df$delta_R2
  # create a name like r2_srage2, r2_lonsmoo, etc.
  root <- sub(".*\\.", "", nm)        # take part after last dot
  col_name <- paste0("r2_", root)
  
  r2_list[[col_name]] <- r2_val
}

# Combine into a single data frame
gam_r2 <- as.data.frame(r2_list)

# Add the region column at the beginning
gam_r2 <- cbind(region = region, gam_r2)

write.csv(gam_r2, row.names = FALSE, 
          file = file.path(data_dir, "from_gam", 
                           paste0(dataset_name, "_GAM_r2_", DV_of_interest, ".csv")))

# =============================================
# Get derivatives from final models        ####
# =============================================

# Compute derivatives of interest for all models
deriv_models <- setdiff(names(final_models), "cov")

# Compute derivatives for remaining models
deriv_list <- setNames(
  lapply(deriv_models, function(name) {
    cat("Computing derivatives for", name, "\n")
    ml <- final_models[[name]]
    
    all_terms_flag <- (name == "sr.mod.srsmoo")
    
    get.derivs(ml, label_of_interest, all_terms = all_terms_flag)
  }),
  deriv_models
)

deriv_df <- dplyr::bind_rows(deriv_list, .id = "model_name")

# Combine by column-binding, keeping only one Parcel column
derivs_all <- cbind(
  region = deriv_list[[1]]$DV,
  do.call(cbind, lapply(deriv_list, function(x) x[, -1]))
)

write.csv(derivs_all, row.names = FALSE, 
          file = file.path(data_dir, "from_gam", 
                           paste0(dataset_name, "_mean_deriv_", DV_of_interest, ".csv")))

