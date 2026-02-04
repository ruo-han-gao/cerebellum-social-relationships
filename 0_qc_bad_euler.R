library(rstudioapi)       # For getting script path
library(dplyr)

# ======================================================
#                 Environment setup                 ####
# ======================================================
rm(list = ls())            # Clear environment
options(scipen = 999, digits = 3) # Disable scientific notation, set digits

# Directories
codedir = dirname(getActiveDocumentContext()$path)
input_dir = "/data/p_03120/hcp_aging_stats/"
output_dir =  "/data/p_03120/hcp_aging_GIN/"

setwd(input_dir)

# ======================================================
#                 List files                        ####
# ======================================================

# Find all aseg.stats files recursively
aseg_files <- list.files(
  path = input_dir,
  pattern = "aseg\\.stats$",    # exact match to aseg.stats
  recursive = TRUE,
  full.names = TRUE
)

# ======================================================
#            Extract aseg info                      ####
# ======================================================

extract_aseg_info <- function(file_path) {
  lines <- readLines(file_path)
  
  # Subject ID (line 5)
  subject <- regmatches(lines[5], regexpr("HCA[0-9]+", lines[5]))
  
  # Surface holes (line 34)
  nholes <- as.numeric(sub(".*, ([0-9]+).*", "\\1", lines[34]))
  
  # ICV (line 35)
  icv <- as.numeric(sub(".*, ([0-9.]+).*", "\\1", lines[35]))
  
  data.frame(
    subject = subject,
    nholes = nholes,
    icv = icv,
    stringsAsFactors = FALSE
  )
}

# Read and combine all files
aseg_all <- do.call(rbind, lapply(aseg_files, extract_aseg_info))

# ICV
icv_df = aseg_all %>% select(-nholes)
write.csv(icv_df, file.path(output_dir, "hcpa_icv.csv"), row.names = FALSE)

# ======================================================
#         Plot distribution of surface holes        ####
# ======================================================
png(file.path(input_dir, "surface_holes_hist_0.png"), width = 1600, height = 1000, res = 300)

hist(
  aseg_all$nholes[aseg_all$nholes],   # Exclude zeros if needed
  main = "Distribution of surface holes",
  xlab = "Total number of defect holes",
  ylab = "Count",
  col = "lightblue",
  border = "black"
)

dev.off()  # Save figure

# ======================================================
#                 Identify outliers
# ======================================================

# try:
m = mean(aseg_all$nholes)
sd = sd(aseg_all$nholes)

m + 2*sd
m + 3*sd

# ===========
med = median((aseg_all$nholes))
mad = mad(aseg_all$nholes)

check = med + 2*mad
check = med + 3*mad
med + 4*mad

##########
hcpa_threshold <- 80   # Based on distribution tail (we defined this post-hoc after plotting the dataset distribution)

euler_bad <- subset(aseg_all, nholes > hcpa_threshold)

# Save subjects with too many surface holes
write.csv(euler_bad, file.path(output_dir, "hcpa_bad_euler_subs.csv"), row.names = FALSE)

# check 
hcpa_threshold <- check   

euler_bad_check <- subset(aseg_all, nholes > hcpa_threshold)

# Save subjects with too many surface holes
write.csv(euler_bad_check, file.path(output_dir, "hcpa_bad_euler_subs_check.csv"), row.names = FALSE)

