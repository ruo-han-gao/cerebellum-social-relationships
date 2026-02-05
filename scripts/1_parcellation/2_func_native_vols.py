#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: manoli and gao
"""
import os
import glob
import numpy as np
import pandas as pd
import nibabel as nib

# Set project directory
project_dir = '/data/p_03120/hcp_aging_GIN'

# Atlas names and corresponding file names
atlas_names = [
    "fusion",
    "MDTB",
    "rest"
]

atlas_files = [
    "fusion/atl-NettekovenSym32_space-MNI_dseg.nii",
    "MDTB/atl-MDTB10_space-MNI_dseg.nii",
    "rest/atl-Buckner7_space-MNI_dseg.nii"
]

# Import subject list
with open(os.path.join(project_dir, 'hcpa_sub_list.txt'), 'r') as file:
    subject_list = [line.strip() for line in file]

# Loop through each atlas
for atlas_name, atlas_file in zip(atlas_names, atlas_files):
    print(f"\nProcessing atlas: {atlas_name}")

    # Load atlas image
    atlas_path = os.path.join(project_dir, "atlases", atlas_file)
    atlas_img = nib.load(atlas_path)
    atlas_data = atlas_img.get_fdata()
    print("Unique labels in atlas:", np.unique(atlas_data))

    # Load region names
    label_file = os.path.join(project_dir, "atlases", f"{atlas_name}", f"{atlas_name}_labels.txt")
    vol_names = pd.read_csv(label_file, header=None)
    vol_names.columns = ['name']
    print(vol_names)

    # Initialize an empty DataFrame to store native space volumes
    vols_df = pd.DataFrame(columns=['subject', 'label', 'volume'])

    # Loop through each subject
    for subject in subject_list:
        mask_dir = os.path.join(project_dir, 'segmentations', atlas_name, subject)
        mask_files = glob.glob(os.path.join(mask_dir, f"{atlas_name}_native.nii*"))

        if not mask_files:
            print(f"Mask not found for subject {subject} in atlas {atlas_name}")
            continue

        # Load mask
        mask_img = nib.load(mask_files[0])
        mask_data = np.round(mask_img.get_fdata(dtype=np.float32)).astype(int)
        voxel_volume = np.prod(mask_img.header.get_zooms())
        unique_labels = np.unique(mask_data)

        parcel_volumes = {
            label: round(np.sum(mask_data == label) * voxel_volume, 2)
            for label in unique_labels if label != 0 # 0 represents background
        }

        vols = pd.DataFrame.from_dict(parcel_volumes, orient='index', columns=['volume'])
        vols.reset_index(inplace=True)
        vols.rename(columns={'index': 'label'}, inplace=True)
        vols['subject'] = subject
        vols = pd.concat([vols, vol_names], axis=1)
        vols_df = pd.concat([vols_df, vols], ignore_index=True)

    # Check for missing values
    missing = vols_df[vols_df['volume'].isnull()]['subject'].unique()
    if len(missing) > 0:
        print(f"Subjects with missing volume data in {atlas_name}:", missing)

    # Convert to wide format
    wide_df = vols_df.pivot_table(index='subject', columns='name', values='volume', aggfunc='first')
    wide_df.reset_index(inplace=True)

    # Save full volume data
    output_dir = os.path.join(project_dir, 'native_space_data')
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"all_{atlas_name}_native_vols.csv")
    wide_df.to_csv(output_file, index=False)
    print(f"Saved volume data: {output_file}")

    # Outlier detection
    mean = wide_df.iloc[:, 1:].mean()
    std_dev = wide_df.iloc[:, 1:].std()
    outlier_threshold_upper = mean + 2 * std_dev
    outlier_threshold_lower = mean - 2 * std_dev

    # Count number of parcels based on the label file
    n_parcels = len(vol_names)
    
    # Count the number of outliers for each subject
    outliers = (wide_df.iloc[:, 1:] > outlier_threshold_upper) | (wide_df.iloc[:, 1:] < outlier_threshold_lower)
    wide_df['outlier_count'] = outliers.sum(axis=1)
    
    # Define outlier threshold per subject (≥ 50% of parcels)
    outlier_cutoff = n_parcels // 2 
    
    # Identify subjects with excessive outliers
    flagged_subjects = wide_df[wide_df['outlier_count'] >= outlier_cutoff]
    print(f"{atlas_name}: {n_parcels} parcels, flagging subjects with ≥ {outlier_cutoff} outliers")

    # Save filtered dataset
    filtered_df = wide_df[~wide_df['subject'].isin(flagged_subjects['subject'])]
    filtered_df = filtered_df.drop(columns=['outlier_count'])
    filtered_df.reset_index(drop=True, inplace=True)

    filtered_file = os.path.join(output_dir, f"{atlas_name}_native_vols_clean.csv")
    filtered_df.to_csv(filtered_file, index=False)
    print(f"Saved filtered data (no outliers): {filtered_file}")

