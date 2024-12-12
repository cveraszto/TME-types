import os, re
import pandas as pd
import random
import numpy as np
import nrrd


def get_all_filenames(folder_path):
    filenames = []
    for filename in os.listdir(folder_path):
        full_path = os.path.join(folder_path, filename)
        if os.path.isfile(full_path):
            filenames.append(filename)
    return filenames

def get_csv_filenames(folder_path):
    csv_filenames = []
    for filename in os.listdir(folder_path):
        full_path = os.path.join(folder_path, filename)
        if os.path.isfile(full_path) and filename.endswith('.csv'):
            csv_filenames.append(filename)
    return sorted(csv_filenames)

def extract_prefix_from_filenames(csv_filenames):
    prefixes = []
    for filename in csv_filenames:
        prefix = filename.split('_')[0]
        prefixes.append(prefix)
    return prefixes

def read_csv_files_into_single_dataframe(filenames, unique_prefixes, folder_path):
    result_dataframes = {}
    
    for prefix in unique_prefixes:
        matching_files = [filename for filename in filenames if filename.startswith(prefix + '_')]
        if matching_files:
            dfs = [pd.read_csv(os.path.join(folder_path, filename), index_col='cluster') for filename in matching_files]
            result_dataframes[prefix] = pd.concat(dfs, ignore_index=True)
    
    return result_dataframes

def read_and_concat_csv_files(filenames, unique_prefixes, folder_path):
    result_dataframes = {}
    
    for prefix in unique_prefixes:
        matching_files = [filename for filename in filenames if filename.startswith(prefix + '_')]
        if matching_files:
            dfs = []
            for filename in matching_files:
                df = pd.read_csv(os.path.join(folder_path, filename), index_col='cluster')
                
                # Check if the DataFrame has 'average_density' and 'density' columns
                if 'average_density' in df.columns:
                    # Keep only 'cluster', 'average_density', and rename 'average_density' to 'concatenated_density'
                    df = df[['average_density']].rename(columns={'average_density': 'concatenated_density'})
                elif 'density' in df.columns:
                    # Keep only 'cluster' and rename 'density' to 'concatenated_density'
                    df = df[['density']].rename(columns={'density': 'concatenated_density'})
                
                dfs.append(df)
            
            result_dataframes[prefix] = pd.concat(dfs, ignore_index=False)
    
    return result_dataframes


def read_and_concat_csv_files_new(filenames, unique_prefixes, folder_path):
    result_dataframes = {}
    
    for prefix in unique_prefixes:
        matching_files = [filename for filename in filenames if filename.startswith(prefix + '_')]
        if matching_files:
            #print(matching_files)
            dfs = []
            for filename in matching_files:
                # Specify the file path
                file_path = os.path.join(folder_path, filename)
                # Read only the 'density_mm3' column into a DataFrame with 'cluster' as the index
                df = pd.read_csv(file_path, usecols=['cluster', 'density_mm3'], index_col='cluster')
                 
                dfs.append(df)
                #print(len(dfs))
    
            result_dataframes[prefix] = pd.concat(dfs, ignore_index=False)
    return result_dataframes


def combine_rows_and_calculate_average(result_dataframes):
    combined_dataframes = {}
    
    for prefix, df in result_dataframes.items():
        combined_dataframes[prefix] = df.groupby(level=0).mean()
    
    return combined_dataframes


def create_combined_dataframe(result_dataframes):
    combined_dataframes = {}
    
    # Get the set of all clusters present in the result_dataframes
    all_clusters = set()
    for df in result_dataframes.values():
        all_clusters.update(df.index)
    
    # Shuffle through the dataframes and select a random cluster from each dataframe
    for cluster in all_clusters:
        selected_dfs = []
        for prefix, df in result_dataframes.items():
            if cluster in df.index:
                selected_row = df.loc[[cluster]].copy()
                # Rename the row index to the unique prefix item
                selected_row.index = [prefix]
                selected_dfs.append(selected_row)
        
        # Combine selected rows into a new DataFrame
        combined_dataframes[cluster] = pd.concat(selected_dfs)
    
    return combined_dataframes


# def nrrd_from_df(df, description_of_all_indexes, CCFv3_0):
#     all_ids_for_df = []
#     df_comb = pd.DataFrame()

#     for regionname in df.index.values:
#         density = df.loc[regionname][0]
#         annotation_id_info = description_of_all_indexes[description_of_all_indexes['parcellation_term_acronym'] == regionname]
#         Annotation2020ids = [int(re.search(r'\d+$', s).group()) for s in annotation_id_info['parcellation_label'].values]
#         df_sub = pd.DataFrame({'density': density}, index=Annotation2020ids)
#         df_comb = pd.concat([df_comb, df_sub])

#         all_ids_for_df.append(Annotation2020ids)

#         print(regionname, annotation_id_info.shape[0], "(region name and its subregions)", flush=True)

#     all_ids_for_df = [value for sublist in all_ids_for_df for value in sublist]
#     all_ids_for_df.append(0)

#     outside = 0
#     outsideid = [0]
#     df_sub = pd.DataFrame({'density': outside}, index=outsideid)
#     df_comb = pd.concat([df_comb, df_sub])

#     # Create float copy of the annotation voluem, update values in CCFv3_0 based on conditions
#     CCFv3_0_copy = np.copy(CCFv3_0).astype('float64')
#     # Expression is 0 here
#     CCFv3_0_copy[~np.isin(CCFv3_0_copy, all_ids_for_df)] = 0.0 
    
#     # Expression is non-zero in these leaf region(s)
#     for index, row in df_comb.iterrows():
#         density_value = row['density']
#         region_id = index
#         CCFv3_0_copy[np.isin(CCFv3_0, region_id)] = density_value
        
#     #Create outside of the brain as np.nan
#     CCFv3_0_copy[np.isin(CCFv3_0, int(0))] = np.nan
    
#     return CCFv3_0_copy

def nrrd_from_df(df, description_of_all_indexes, CCFv3_0):
    """
    Creates a copy of the CCFv3, 3D numpy array, where the voxel values are of densities of one cell-type.
    Relies on global variable description_of_all_indexes which is an edited version of the original Allen's
    parcellation_to_parcellation_term_membership file. 
    Parameters:
    - df: DataFrame, with one column of densities, where rows are brain regions.
    - description_of_all_indexes: df of information on brain areas from the Allen Institute
    - CCFv3_0: 3D numpy array representing a version of the Allen annotation volume .

    Returns:
    - CCFv3_0_copy: 3D numpy array of the annotation volume where voxels are densities of a cell-type,
    while outside the brain is set to np.nan. 
    """
   
    all_ids_for_df = []
    df_comb = pd.DataFrame()

    for regionname in df.index.values:
        density = df.loc[regionname][0]
        annotation_id_info = description_of_all_indexes[description_of_all_indexes['cluster_as_filename'] == regionname]
        Annotation2020ids = list(annotation_id_info['label_numbers'].values)
        df_sub = pd.DataFrame({'density': density}, index=Annotation2020ids)
        df_comb = pd.concat([df_comb, df_sub])

        all_ids_for_df.append(Annotation2020ids)

        if annotation_id_info.shape[0] == 0:
            print(regionname, " region was not found and the its density will be 0!!!", flush=True)
        elif annotation_id_info.shape[0] == 1:
            print(regionname, annotation_id_info.shape[0], "(region name is a leaf region )", flush=True)
        else:
            print(regionname, annotation_id_info.shape[0], "(region name and its subregions)", flush=True)

    #After the loop: only add the outside 0 once
    all_ids_for_df = [value for sublist in all_ids_for_df for value in sublist]
    all_ids_for_df.append(0)

    outside = 0
    outsideid = [0]
    df_sub = pd.DataFrame({'density': outside}, index=outsideid)
    df_comb = pd.concat([df_comb, df_sub])

    # Create float copy of the annotation volume, update values in CCFv3_0 based on conditions
    CCFv3_0_copy = np.copy(CCFv3_0).astype('float64')
    # Expression is 0 here
    CCFv3_0_copy[~np.isin(CCFv3_0_copy, all_ids_for_df)] = 0.0 
    
    # Expression is non-zero in these leaf region(s)
    for index, row in df_comb.iterrows():
        density_value = row['density']
        region_id = index
        CCFv3_0_copy[np.isin(CCFv3_0, region_id)] = density_value
        
    #Create outside of the brain as np.nan
    CCFv3_0_copy[np.isin(CCFv3_0, int(0))] = np.nan
    
    return CCFv3_0_copy

def get_nrrd_files(folder_path):
    """
    Get NRRD files from a folder and return the file names with and without extensions.

    Parameters:
    - folder_path: Path to the folder containing NRRD files.

    Returns:
    - nrrd_files: List of NRRD file names.
    - nrrd_files_without_extension: List of NRRD file names without extensions.
    """
    # Get all files in the folder
    files = os.listdir(folder_path)

    # Filter files with the .nrrd extension
    nrrd_files = [file for file in files if file.endswith('.nrrd')]

    # Get NRRD file names without extension
    nrrd_files_without_extension = [file_name[:-5] for file_name in nrrd_files]

    return nrrd_files, nrrd_files_without_extension

#Load all csv files with densities
root_folder = '/gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/results/density_calculations/'
folder_path = f'{root_folder}csv/'
filenames = get_all_filenames(folder_path)
csv_filenames = get_csv_filenames(folder_path)
prefixes = extract_prefix_from_filenames(csv_filenames)
unique_prefixes = sorted(list(set(prefixes)))

#Create a dict of df, each containing a cell type's occurence in all regions and its densities in all regions
result_dataframes = read_and_concat_csv_files_new(csv_filenames, unique_prefixes, folder_path)
combined_result_dataframes = combine_rows_and_calculate_average(result_dataframes)
shuffled_combined_dataframes = create_combined_dataframe(combined_result_dataframes)

# Get the set of all clusters present in the result_dataframes
all_clusters = set()
for df in result_dataframes.values():
    all_clusters.update(df.index)

cluster_ids = []
for cluster in all_clusters: 
    match = re.search(r'\d+', cluster)
    if match:
        # Extract the first 4 numbers
        cluster = match.group()[:4]
    else:
        # If no numeric portion is found, keep the original string
        cluster = cluster
        print(f"No numeric portion is found, keeping original name for: {cluster}", flush=True)
    
    cluster_ids.append(cluster)
     
cluster_ids = sorted(cluster_ids)

#region_list = list(combined_result_dataframes.keys())

#This block is to get all region names for placement into CCFv3 space
view_directory = '/gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/metadata/MERFISH-C57BL6J-638850-CCF/20231215/views'
file = os.path.join( view_directory, 'cell_metadata_with_parcellation_annotation.csv')
cell_joined = pd.read_csv(file)
cell_joined.set_index('cell_label',inplace=True)
#file = '/gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/metadata/Allen-CCF-2020/20230630/parcellation_to_parcellation_term_membership.csv'
file = '/gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/metadata/parcellation_to_parcellation_term_membership_extend.csv'
parcellation_annotation = pd.read_csv(file)
#This describes all indexes in the slices, but also on every level: organ category division structure substructure
parcellation_indexes = list(np.unique(cell_joined['parcellation_index']))
description_of_all_indexes = parcellation_annotation[parcellation_annotation['parcellation_index'].isin(parcellation_indexes)]

#Read CCFv3 annotation volumes (choose 1)
# Specify the full or relative path to the log file
# #1
# data_folder = "/gpfs/bbp.cscs.ch/project/proj84/piluso/share/general/warped_augmented_CCFv3/"
# CCFv3_0, _ = nrrd.read(f'{data_folder}annotation_25_2022_CCFv3_0.nrrd')
# save_nrrd = f'{root_folder}nrrd_CCFv3_0/'
# #2
# data_folder = "/gpfs/bbp.cscs.ch/data/project/proj62/csaba/atlas/bbp_prod_files/2022/"
# CCFv3_0, _ = nrrd.read(f'{data_folder}annotation_25.nrrd')
# save_nrrd = f'{root_folder}nrrd/'
#3
data_folder = "/gpfs/bbp.cscs.ch/project/proj84/piluso/share/general/warped_augmented_CCFv3/"
CCFv3_0, _ = nrrd.read(f'{data_folder}annotation_25_2022_CCFv3a.nrrd')
save_nrrd = f'{root_folder}nrrd_CCFv3a/'
# #4 This is the BB version
# CCFv3_0, _ = nrrd.read("/gpfs/bbp.cscs.ch/data/project/proj84/atlas_pipeline_runs/2024-05-15T22:44:26+02:00/annotation_ccfv3_l23split_barrelsplit_validated.nrrd")
# save_nrrd = f'{root_folder}scaled_nrrd_CCFv3a_validated/'

#Main loop running on each node
count = 0
for cluster, df in shuffled_combined_dataframes.items():
    print(f"Combined DataFrame for cluster '{cluster}':", flush=True)
    f_name = cluster.translate(str.maketrans({"/":  r"", "-":  r"_", " ":  r"_"}))
    created_nrrds, nrrd_files_without_extension = get_nrrd_files(save_nrrd)
    if f_name not in nrrd_files_without_extension:
        print(f"{f_name} is not created yet.", flush=True)

        #code
        result_CCFv3_0_copy = nrrd_from_df(df, description_of_all_indexes, CCFv3_0)

        print(f"Saving file: {cluster} as {f_name}" )
        print("\n", flush=True)
        nrrd.write(f"{save_nrrd}{f_name}.nrrd", result_CCFv3_0_copy)
        #iterate      
        count += 1
        print(count, flush=True)
        # if count == 1:
        #     break

    else:
        print(f"Exception occurred for cluster {cluster}. File already exists.", flush=True)
        # Handle the exception if needed
        print("\n", flush=True)
        continue