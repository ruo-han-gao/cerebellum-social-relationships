library(rstudioapi) # gets path of script
library(dplyr)
library(stringr)
library(purrr)

# clear the environment
rm(list = ls())
options(scipen = 999, digits = 3)

codedir = dirname(getActiveDocumentContext()$path)
outputdir = "/data/tu_gaor/scripts/dat"

hcpd_dir = "/data/p_03120/"
hcpa_dir = "/data/p_03120/hcp_aging_GIN/"

# specify source data dir here
source_dir = hcpd_dir

if (source_dir == hcpa_dir) {
  dataset_name <- "hcpa"
} else if (source_dir == hcpd_dir) {
  dataset_name <- "hcpd"
} else {
  dataset_name <- "unknown"
}

# build the dataset-specific output dir
output_dir <- file.path(outputdir, dataset_name)

setwd(source_dir)

# =============================================
#   Rename hcpd file the same as hcpa      ####
# =============================================

if (source_dir == hcpd_dir) {
  bad_euler = read.csv(file.path(source_dir, "data/bad_euler_subs_hcpd.csv"))
  bad_euler$subject = sub("_V1_MR$", "", bad_euler$subject_ID)
  write.csv(bad_euler, file.path(source_dir, "hcpd_bad_euler_subs.csv"))
  discarded = read.table(file.path(source_dir, "data/hcpd_discarded_subs.txt"))
  discarded$subject = sub("_V1_MR$", "", discarded$V1)
  write.csv(discarded, file.path(source_dir, "hcpd_discarded_subs.csv"))
} else {
  message("Skipping renaming files for hcpa")
}

# =============================================
#        Native cerebellum data            ####
# =============================================

fusion = read.csv(file.path(source_dir, "native_space_data/fusion_native_vols_clean.csv"))
# mdtb = read.csv(file = "native_space_data/MDTB_native_vols_clean.csv")
# rest = read.csv(file = "native_space_data/rest_native_vols_clean.csv")

cerebel_data <- fusion


# =============================================
#        Native cortical data            ####
# =============================================

LH_yeo7 = read.delim(file.path(source_dir, paste0("data/", dataset_name, "_Yeo7_lh_volume.tsv")))
RH_yeo7 = read.delim(file.path(source_dir, paste0("data/", dataset_name, "_Yeo7_rh_volume.tsv")))

LH_yeo17 = read.delim(file.path(source_dir, paste0("data/", dataset_name, "_Yeo17_lh_volume.tsv")))
RH_yeo17 = read.delim(file.path(source_dir, paste0("data/", dataset_name, "_Yeo17_rh_volume.tsv")))

rename_yeo <- function(df) {
  
  # detect hemisphere prefix from the first column
  first_col <- colnames(df)[1]
  hemi <- if (grepl("^rh\\.", first_col)) {
    "rh"
  } else if (grepl("^lh\\.", first_col)) {
    "lh"
  } else {
    stop("Could not detect hemisphere prefix (expected rh. or lh.)")
  }
  
  # build regex patterns
  subject_pattern <- paste0("^", hemi, "\\.Yeo2011_[0-9]+Networks_N1000\\.volume$")
  network_pattern <- paste0("^", hemi, "_[0-9]+Networks_[0-9]+_volume$")
  
  # detect actual columns
  subject_col <- grep(subject_pattern, colnames(df), value = TRUE)
  network_cols <- grep(network_pattern, colnames(df), value = TRUE)
  
  # subset
  df_sub <- df[, c(subject_col, network_cols)]
  
  # rename columns
  new_names <- c("subject", gsub("_volume$", "", network_cols))
  colnames(df_sub) <- new_names
  
  # clean subject ID
  if (any(grepl("_V1_MR$", df_sub$subject))) {
    df_sub$subject <- sub("_V1_MR$", "", df_sub$subject)
  }
  
  return(df_sub)
}

yeo7_lh = rename_yeo(LH_yeo7)
yeo7_rh = rename_yeo(RH_yeo7)
write.csv(yeo7_lh, file.path(source_dir, paste0(dataset_name, "_Yeo7_lh_volume.csv")))
write.csv(yeo7_rh, file.path(source_dir, paste0(dataset_name, "_Yeo7_rh_volume.csv")))

yeo17_lh = rename_yeo(LH_yeo17)
yeo17_rh = rename_yeo(RH_yeo17)
write.csv(yeo17_lh, file.path(source_dir, paste0(dataset_name, "_Yeo17_lh_volume.csv")))
write.csv(yeo17_rh, file.path(source_dir, paste0(dataset_name, "_Yeo17_rh_volume.csv")))

yeo7  = yeo7_lh  %>% left_join(yeo7_rh,  by = "subject")
yeo17 = yeo17_lh %>% left_join(yeo17_rh, by = "subject")

reorder_hemis <- function(df, hemi_prefixes = c("lh", "rh")) {
  # keep subject first
  subject_col <- "subject"
  
  # detect left and right columns
  lh_cols <- grep("^lh_", colnames(df), value = TRUE)
  rh_cols <- grep("^rh_", colnames(df), value = TRUE)
  
  # extract network numbers to sort properly
  get_num <- function(x) as.numeric(sub(".*_(\\d+)$", "\\1", x))
  lh_cols <- lh_cols[order(get_num(lh_cols))]
  rh_cols <- rh_cols[order(get_num(rh_cols))]
  
  # interleave lh and rh
  interleaved <- as.vector(rbind(lh_cols, rh_cols))
  
  # final column order
  df <- df[, c(subject_col, interleaved)]
  
  return(df)
}

yeo7_reordered <- reorder_hemis(yeo7)
yeo17_reordered <- reorder_hemis(yeo17)

# =============================================
#             Prep functions               ####
# =============================================

sum_data_available <- function(df) {
  n_subjects <- length(unique(df$src_subject_id))
  age_range <- range(df$interview_age, na.rm = TRUE)/12
  
  number_of_values <- sapply(df, function(x) sum(!is.na(x) & x != ""))
  
  cat("Number of subjects:", n_subjects, "\n")
  cat(sprintf("Age range: %.2f to %.2f years\n", age_range[1], age_range[2]), "\n")
  cat("Number of values per column:\n")
  print(number_of_values)
}

sum_items <- function(df, selected_cols = NULL) {
  df_sub <- df[, selected_cols, drop = FALSE]
  values <- trimws(as.vector(unlist(df_sub[1, ])))
  item_df <- data.frame(
    rowname = colnames(df_sub),
    value = values,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  return(item_df)
}

# =============================================
#           Emotional support              ####
# =============================================

emosup = read.table("phenotype/tlbx_emsup01.txt", header = TRUE)

emosup_scores <- emosup %>%
  select(src_subject_id, interview_age, 
         soc200, soc203, soc205, soc216, soc222)%>%
  slice(-1) %>% 
  mutate(across(-src_subject_id, as.numeric))

emosup_scores0 <- emosup_scores %>%
  group_by(src_subject_id) %>%
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop") %>%
  mutate(raw_emosup = rowSums(across(3:last_col())))

sum_data_available(emosup_scores0)
emo_both = c("soc200", "soc203", "soc205","soc216", "soc222")
emo_items = sum_items(emosup, emo_both)

# =============================================
#                Friendship                ####
# =============================================

friends = read.table("phenotype/tlbx_friend01.txt", header = TRUE)
friends_scores <- friends %>%
  select(src_subject_id, interview_age, 
         soc230, soc233, soc237, soc239, soc247) %>%
  slice(-1) %>% 
  mutate(across(-src_subject_id, as.numeric)) %>%
  mutate(raw_friends = rowSums(across(3:last_col())))

sum_data_available(friends_scores)

friends_both = c("soc230", "soc233", "soc237", "soc239", "soc247")
friends_items = sum_items(friends, friends_both)

# =============================================
#                Loneliness                ####
# =============================================

loneli = read.table("phenotype/prsi01.txt", header = TRUE, sep = "\t")
loneli = loneli %>%
  rename(soc254 = ucla11x2)

loneli_scores <- loneli %>%
  select(src_subject_id, interview_age, 
         soc253, soc254, soc260, soc261) %>%
  slice(-1) %>% 
  mutate(across(-src_subject_id, as.numeric)) %>%
  mutate(raw_loneli = rowSums(across(3:last_col())))

sum_data_available(loneli_scores)

loneli_both = c("soc253", "soc254", "soc260", "soc261")
lonely_items = sum_items(loneli, loneli_both)

# =============================================
#           Perceived rejection            ####
# =============================================

rej = read.table("phenotype/tlbx_rej01.txt", header = TRUE)
# table(rej$version_form)

per_rej_scores = rej %>%
  filter(collection_id == "collection_id" | 
           rej$version_form == "NIH Toolbox Perceived Rejection FF Age 18+ v2.0" | 
           rej$version_form == "NIH Toolbox Perceived Rejection FF Ages 8-17 v2.0") %>%
  select(src_subject_id, interview_age, 
         soc276, soc279, soc281) %>%
  slice(-1) %>% 
  mutate(across(-src_subject_id, as.numeric)) %>%
  mutate(raw_rej = rowSums(across(3:last_col())))

rej_both = c("soc276", "soc279", "soc281")
sum_data_available(per_rej_scores)
rej_items = sum_items(rej, rej_both)

# =============================================
#           Perceived hostility            ####
# =============================================

hostility = read.table("phenotype/tlbx_perhost01.txt", header = TRUE, sep = "\t")

hostility_scores = hostility %>%
  select(src_subject_id, interview_age,
         soc262, soc263, soc267, soc268, soc270)%>%
  slice(-1) %>% 
  mutate(across(-src_subject_id, as.numeric))

hostility_both = c("soc262", "soc263", "soc267", "soc268", "soc270")
sum_data_available(hostility_scores)
hostility_items = sum_items(hostility, hostility_both)


#impute missing values
if (source_dir == hcpd_dir) {
  rows_missing <- which(rowSums(is.na(hostility_scores)) > 0)
  cols_to_impute <- 3:7
  
  for (row_num in rows_missing) {
    # Compute median of non-missing values in the selected columns
    med_for_missing <- median(unlist(hostility_scores[row_num, cols_to_impute]), na.rm = TRUE)
    
    # Replace NAs with the median
    hostility_scores[row_num, cols_to_impute][is.na(hostility_scores[row_num, cols_to_impute])] <- med_for_missing
  }
} else {
  message("Skipping missing value imputation for hcpa")
}
hostility_scores = hostility_scores %>%
  mutate(raw_hostility = rowSums(across(3:last_col())))
sum_data_available(hostility_scores)

# =============================================
#            Combine social scales         ####
# =============================================

itemname_map <- c(
  emo200 = "soc200", emo203 = "soc203", emo205 = "soc205", emo222 = "soc222", emo216 = "soc216",
  fri230 = "soc230", fri233 = "soc233", fri237 = "soc237", fri239 = "soc239", fri247 = "soc247",
  lon253 = "soc253", lon254 = "soc254", lon260 = "soc260", lon261 = "soc261",
  rej276 = "soc276", rej279 = "soc279", rej281 = "soc281",
  hos262 = "soc262", hos263 = "soc263", hos267 = "soc267", hos268 = "soc268", hos270 = "soc270"
)

# all items ====
combined_items <- bind_rows(
  emo = emo_items,
  friend = friends_items,
  lonely = lonely_items,
  rejection = rej_items,
  hostility = hostility_items,
  .id = "scale"
)
combined_items = combined_items %>% rename(itemnames = rowname)
combined_items <- combined_items %>%
  mutate(itemnames = recode(itemnames, !!!setNames(names(itemname_map), itemname_map)))

# write item content in a file
output_file <- file.path(outputdir, "combined_items.csv")
write.csv(combined_items, output_file, row.names = FALSE)

# all scores ====
df_list = list(friends_scores, emosup_scores0, loneli_scores, per_rej_scores, hostility_scores) %>%
  map(~ select(.x, -interview_age))
social_itemscores <- reduce(df_list, left_join, by = "src_subject_id")
social_itemscores = social_itemscores %>% rename(subject = src_subject_id)

# rename
social_itemscores <- social_itemscores %>%
  rename(!!!itemname_map)

# put raw scores at the end of the df
raw_cols <- c("raw_friends", "raw_emosup", "raw_loneli", "raw_rej", "raw_hostility")
social_itemscores <- social_itemscores %>%
  select(-all_of(raw_cols), all_of(raw_cols))

social_itemscores[which(rowSums(is.na(social_itemscores)) > 0), ]
social_itemscores[, 2:23] <- lapply(social_itemscores[, 2:23], as.numeric)
colMeans(social_itemscores[, 2:28])

range(social_itemscores[,2:23])

# =============================================
#             Valence conversion           ####
# =============================================

cols_for_conver <- social_itemscores %>%
  select(starts_with("lon"), starts_with("rej"), starts_with("hos")) %>%
  colnames()

social_itemscores <- social_itemscores %>%
  mutate(across(all_of(cols_for_conver), ~6 - .))
colMeans(social_itemscores[, 2:28])
range(social_itemscores[,2:23])

# =============================================
#                 Age and sex              ####
# =============================================

# age and sex were taken from the friends scale because for both hcpd and hcpa, this file contains all the subjects and no duplicates
age_sex <- friends %>%
  slice(-1) %>%  # remove the description row
  select(src_subject_id, interview_age, sex) %>%
  rename(subject = src_subject_id) %>%
  mutate(age = as.numeric(interview_age) / 12) %>%
  select(-interview_age)

# =============================================
#        Discard subjects in hcpA           ####
# =============================================

if (source_dir == hcpa_dir) {
  # age == 100
  age100 = age_sex %>%
    filter(age == 100)
  
  # failed in yeo parcellation
  failed = read.table(file.path(source_dir, paste0("data/", dataset_name, "_failed_subjects.txt")))
  colnames(failed) <- "subject"
  if (any(grepl("_V1_MR$", failed$subject))) {
    failed$subject <- sub("_V1_MR$", "", failed$subject)
  }
  discarded = bind_rows(age100, failed)
  write.csv(discarded, file.path(source_dir, "hcpa_discarded_subs.csv"))
} else {
  message("Skipping *discarded file* for hcpd, alredy there")
}

# =============================================
#                     ICV                  ####
# =============================================
icv = read.csv(file.path(source_dir, paste0(dataset_name, "_icv.csv")))


# =============================================
#                 Exclusion                ####
# =============================================
bad_euler = read.csv(file.path(source_dir, paste0(dataset_name, "_bad_euler_subs.csv")))
discarded = read.csv(file.path(source_dir, paste0(dataset_name, "_discarded_subs.csv")))

# =============================================
#          Combine and filter data         ####
# =============================================

# combine participants with all data available 
data <- age_sex %>%
  inner_join(icv, by = "subject") %>%
  inner_join(cerebel_data, by = "subject") %>%
  inner_join(yeo7_reordered, by = "subject") %>%
  inner_join(yeo17_reordered, by = "subject") %>%
  inner_join(social_itemscores, by = "subject")
dim(data)

data <- data %>%
  filter(!subject %in% bad_euler$subject) %>%
  filter(!subject %in% discarded$subject)
dim(data)

nrow_data <- nrow(data)

# write cerebellar fusion data and item scores in a file
output_file <- file.path(output_dir, paste0("fusion_yeo7_17_socialscore_n", nrow_data, ".csv"))
write.csv(data, output_file, row.names = FALSE)



