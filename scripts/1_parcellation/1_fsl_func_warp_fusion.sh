#!/bin/bash

# Set base directories
project_dir="/data/p_03120/hcp_aging_GIN"
input_dir="/data/p_03120/hcp_aging_GIN"
output_dir="${project_dir}/segmentations"

# Subject list
subject_list="${project_dir}/hcpa_sub_list.txt"

# Atlas names 
atlas_names=(
  "fusion"
)

# Atlas files (in same order as atlas names)
atlases=(
  "fusion/atl-NettekovenSym32_space-MNI_dseg.nii"
)

# Loop over each subject
while IFS= read -r subject; do
  echo "Processing subject: ${subject}"
  
  # Loop over atlases and their corresponding output names
  for i in "${!atlases[@]}"; do
    atlas_file="${atlases[$i]}"
    atlas_name="${atlas_names[$i]}"
    
    mkdir -p "${output_dir}/${atlas_name}/${subject}"
    
    # Warp MNI-aligned atlas to native space 
    applywarp \
      --ref="${input_dir}/${subject}/T1w_acpc_dc_restore.nii.gz" \
      --in="${project_dir}/atlases/${atlas_file}" \
      --out="${output_dir}/${atlas_name}/${subject}/${atlas_name}_native.nii.gz" \
      --warp="${project_dir}/${subject}/standard2acpc_dc.nii.gz" \
      --interp=nn
  done

done < "$subject_list"

