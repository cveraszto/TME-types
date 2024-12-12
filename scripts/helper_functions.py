'''
Functions related to density estimation from MERFISH slices 
/gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/notebooks
'''
# helper functions . py

import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
import pandas as pd

def calculate_density_mm3(section, volume_single_cube_mm3, selected_region):
    """
    Calculate density (count per mm3_vol) for each cluster in the given DataFrame section.

    Parameters:
        section (DataFrame): DataFrame containing 'cluster' and 'voxel_vol' columns.
        volume_single_cube_mm3 (float): Volume of a single voxel in cubic millimeters.

    Returns:
        DataFrame: DataFrame containing 'cluster' and 'density_mm3' columns.
    """
    # Calculate density (count per voxel_vol) for each cluster
    data = []
    for cluster, count in section['cluster'].value_counts().items():
        #Extract all slices voxel count once and take the sum since cell count is summed too
        voxel_vol = section[['brain_section_label', 'voxel_vol']].drop_duplicates()['voxel_vol'].sum()
        count_per_mm3_vol = count / (voxel_vol * volume_single_cube_mm3)
        data.append({'cluster': cluster, 'density_mm3': count_per_mm3_vol})

    # Create DataFrame from the calculated data
    result_df = pd.DataFrame(data)
    result_df.index.name = selected_region #Store real cell type name in every df
    #result_df.set_index('cluster', inplace=True)
    return result_df
    

def calculate_intersection_vol(CCFv3_0, region_id, z_index, intersection_angle_deg):
    
    """
    Calculate intersection between a 3D array and a 2D orthogonal plane.
    
    Parameters:
    - CCFv3: 3D numpy array representing a version of the Allen annotation volume .
    - region_id: int, representing the region id, of a brain region in the annotation volume.
    - z_index: int, the coronal index where the the 2D plane first intersects the 3D array
    - intersection_angle_deg: float, one can define an angle different from perpendicular
    Returns:
    - intersection_vol_in_voxels: int, the number of the brain region's voxels which were crossed by the 2D plane, i.e. a volume
    """
    
    CCFv3_0_copy = np.copy(CCFv3_0).astype(float)

    # Convert angle to radians
    intersection_angle_rad = np.radians(intersection_angle_deg)

    # Define dimensions of the volume
    volume_shape = CCFv3_0_copy.shape
    array2D_width, array2D_height = volume_shape[2], volume_shape[1]

    # Create meshgrid for x and y coordinates
    x_coordinates, y_coordinates = np.meshgrid(np.arange(array2D_width), np.arange(array2D_height))

    # Calculate y coordinate of the line at each x coordinate based on intersection angle
    y_line = np.tan(intersection_angle_rad) * np.arange(array2D_width)
    y_line_int_clipped = np.clip(np.round(y_line).astype(int), 0, array2D_height - 1)

    # Create boolean mask to identify voxels in intersection region
    mask = (y_coordinates <= y_line_int_clipped)

    # Multiply voxels in intersection region by 2
    CCFv3_0_copy[z_index, :, :][~mask] = 0
    CCFv3_0_copy[z_index, :, :][mask] *= 2

    # Count number of voxels affected by the plane
    intersection_vol_in_voxels = np.count_nonzero(CCFv3_0_copy[z_index, :, :] == region_id*2)

    return intersection_vol_in_voxels


def calculate_total_area(tri, points):
    """
    Calculate the total area of a Delaunay triangulation.

    Parameters:
    - tri: Delaunay triangulation object
    - points: 2D array, input points for triangulation

    Returns:
    - total_area: float, total area of the Delaunay triangulation
    """
    total_area = 0.0

    for simplex in tri.simplices:
        # Extract the vertices of the triangle
        triangle_vertices = points[simplex]

        # Calculate two vectors representing sides of the triangle
        AB = triangle_vertices[1] - triangle_vertices[0]
        AC = triangle_vertices[2] - triangle_vertices[0]

        # Calculate the cross product of AB and AC
        cross_product = np.cross(AB, AC)

        # Calculate the area of the triangle (half the magnitude of the cross product)
        triangle_area = 0.5 * np.linalg.norm(cross_product)

        # Sum up the areas
        total_area += triangle_area

    return total_area


def combine_rows_and_calculate_average(result_dataframes):
    """
    Take a dict of dataframes and calculates the mean of rows with the same index. 
    Since the input is a combination of multiple slices of the same prefixes, ie. brain regions, 
    we have multiple densities of the same cell types in the dfs, thus the resulting dfs may
    become smaller.

    Parameters:
    - result_dataframe: dict, multiple dataframes of cell-type densities in brain regions

    Returns:
    - combined_dataframes: dict, of dfs of average cell type densities
    """
   
    combined_dataframes = {}
    
    for prefix, df in result_dataframes.items():
        combined_dataframes[prefix] = df.groupby(level=0).mean()
    
    return combined_dataframes


def create_combined_dataframe(result_dataframes):
    """
    Takes a dict of dfs of different cell-type densities in brain regions (prefix) and 
    creates a dict of dfs of cell-type densities (keys) in different brain regions.
    Cell-types are also called clusters.  

    Parameters:
    - result_dataframes: dict, multiple dfs of cell-type densities in brain regions

    Returns:
    - combined_dataframes: dict of dfs where clusters are keys and brain region wise
    densities are the values. 
    """
   
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


def delaunay_triangulation_with_perturbation(points, perturbation=1e-6):
    """
    Perform Delaunay triangulation with perturbation to avoid coplanarity issues.

    Parameters:
    - points: 2D array of shape (n, 3) representing the coordinates of points.
    - perturbation: Small value added to each coordinate to avoid coplanarity.

    Returns:
    - perturbed_points: 2D array of shape (n, 3) representing the perturbed points.
    """

    # Copy the original points to avoid modifying the input array
    perturbed_points = points.copy()

    # Add a small perturbation to avoid coplanarity
    perturbed_points += np.random.uniform(-perturbation, perturbation, size=perturbed_points.shape)

    return perturbed_points

def extract_regions_from_column_names(folder_path, file_list):
    '''Reads filenames in a folder, loads them as df to extract the column name. Removes repetitions.'''
    allen_regions = []
    for filename in file_list:
        file_path = os.path.join(folder_path, filename)
        df = pd.read_csv(file_path, nrows=1)  # Read only the first row to get column names
        allen_regions.append(df.columns[0])
        
    unique_allen_regions = sorted(list(set(allen_regions)))    
    return unique_allen_regions

def extract_prefix_from_filenames(csv_filenames):
    '''Takes the first part of a list of filenames before the first _. Prefix = brain region'''
    prefixes = []
    for filename in csv_filenames:
        prefix = filename.split('_')[0]
        prefixes.append(prefix)
    return prefixes


def get_all_filenames(folder_path):
    '''Get all file names from a folder.'''
    filenames = []
    for filename in os.listdir(folder_path):
        full_path = os.path.join(folder_path, filename)
        if os.path.isfile(full_path):
            filenames.append(filename)
    return filenames


def get_csv_filenames(folder_path):
    '''Filters out .csv files from a list of filenames.'''
    csv_filenames = []
    for filename in os.listdir(folder_path):
        full_path = os.path.join(folder_path, filename)
        if os.path.isfile(full_path) and filename.endswith('.csv'):
            csv_filenames.append(filename)
    return sorted(csv_filenames)

def get_csv_files(folder_path):
    """
    Get _ files from a folder and return the file names with and without extensions.

    Parameters:
    - folder_path: Path to the folder containing _ files.

    Returns:
    - cs_files: List of _ file names.
    - csv_files_without_extension: List of _ file names without extensions.
    """
    # Get all files in the folder
    files = os.listdir(folder_path)

    # Filter files with the .nrrd extension
    csv_files = [file for file in files if file.endswith('.csv')]

    # Get CSV file names without extension
    csv_files_without_extension = [file_name[:-5] for file_name in csv_files]

    return csv_files, csv_files_without_extension

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


def nrrd_from_df(df, description_of_all_indexes, CCFv3_0, region_map):
    """
    Creates a copy of the CCFv3, 3D numpy array, where the voxel values are of densities of one cell-type.
    Relies on global variable description_of_all_indexes which is an edited version of the original Allen's
    parcellation_to_parcellation_term_membership file. 
    Parameters:
    - df: DataFrame, with one column of densities, where rows are brain regions.
    - description_of_all_indexes: df of information on brain areas from the Allen Institute
    - CCFv3_0: 3D numpy array representing a version of the Allen annotation volume.
    - region_map: annotation hierarchy loaded from json file (from AIBS)

    Returns:
    - CCFv3_0_copy: 3D numpy array of the annotation volume where voxels are densities of a cell-type,
    while outside the brain is set to np.nan. 
    """
   
    all_ids_for_df = []
    df_comb = pd.DataFrame()

    for regionname in df.index.values:
        density = df.loc[regionname, 'density_mm3']
        annotation_id_info = description_of_all_indexes[description_of_all_indexes['cluster_as_filename'] == regionname]
        Annotation2020ids = list(annotation_id_info['label_numbers'].values) + list(region_map.find(annotation_id_info['parcellation_term_name'].values[0], "name", with_descendants=True))
        df_sub = pd.DataFrame({'density': density}, index=Annotation2020ids)
        df_comb = pd.concat([df_comb, df_sub])

        all_ids_for_df.append(Annotation2020ids) #This is a list only to see where will be cell types in the nrrd

        if annotation_id_info.shape[0] == 0:
            print(regionname, " region was not found and the its density will be 0!!!", flush=True)
        elif annotation_id_info.shape[0] == 1:
            print(regionname, annotation_id_info.shape[0], "(region name is a leaf region )", flush=True)
        else:
            print(regionname, annotation_id_info.shape[0], "(region name and its subregions)", flush=True)

    #After the loop: only add the outside 0 once
    all_ids_for_df = [value for sublist in all_ids_for_df for value in sublist]
    all_ids_for_df.append(0)
    all_ids_for_df = list(set(all_ids_for_df)) #Remove duplicates

    outside = 0
    outsideid = [0]
    df_sub = pd.DataFrame({'density': outside}, index=outsideid)
    df_comb = pd.concat([df_comb, df_sub])
    #Drop Duplicate rows coming from "parcellation duplication"
    df_comb = df_comb.drop_duplicates()

    # Create float copy of the annotation volume, update values in CCFv3_0 based on conditions
    CCFv3_0_copy = np.copy(CCFv3_0).astype('float64')
    # Set expression to 0 (where the cell types don't exist)
    CCFv3_0_copy[~np.isin(CCFv3_0_copy, all_ids_for_df)] = 0.0 
    
    # Expression is non-zero in these leaf region(s)
    for index, row in df_comb.iterrows():
        density_value = row['density']
        region_id = index
        CCFv3_0_copy[np.isin(CCFv3_0, region_id)] = density_value
        
    #Create outside of the brain as np.nan
    CCFv3_0_copy[np.isin(CCFv3_0, int(0))] = np.nan
    
    return CCFv3_0_copy


def plot_section(xx=None, yy=None, ax=None, cc=None, val=None, pcmap=None, 
                 overlay=None, extent=None, bcmap=plt.cm.Greys_r, alpha=1.0,
                 fig_width=6, fig_height=6):
    
    if ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(fig_width, fig_height)

    if xx is not None and yy is not None and pcmap is not None:
        ax.scatter(xx, yy, s=20, c=val, marker='.', cmap=pcmap)
    elif xx is not None and yy is not None and cc is not None:
        ax.scatter(xx, yy, s=20, color=cc, marker='.', zorder=1)
        
    if overlay is not None and extent is not None and bcmap is not None:
        ax.imshow(overlay, cmap=bcmap, extent=extent, alpha=alpha, zorder=2)
        
    ax.set_ylim(11, 0)
    ax.set_xlim(0, 11)
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    
    return ax

def plot_and_save_delaunay(points, tri, title, save_path):
    # Plot the Delaunay triangulation
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_trisurf(points[:, 0], points[:, 1], points[:, 2], triangles=tri.simplices, alpha=0.2)

    # Optionally, you can also plot the original points
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='r', marker='o')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Set the title
    fig.suptitle(title)

    # Save the image instead of showing it
    plt.savefig(save_path)
    plt.close()  # Close the figure to free up resources

    
def process_slide(section, select_section, selected_region, download_base):
    print(f"Number of cells on one-sided slide {select_section} of {selected_region} region is {section.shape[0]}.", file=log_file)          

    # Convert DataFrame to NumPy array
    points = section[['x_reconstructed', 'y_reconstructed', 'z_reconstructed']].values

    # Perform Delaunay triangulation and save image
    perturbed_points = delaunay_triangulation_with_perturbation(points)
    tri = Delaunay(perturbed_points)
    title = f"{selected_region}_SLIDE{np.round(select_section, 2)}_one-sided"
    #cluster_as_folder = selected_region.translate(str.maketrans({"/":  r"", "-":  r"", " ":  r"", ",": r""}))
    save_path = f"{download_base}results/density_calculations_1/img/{cluster_as_folder}_{np.round(select_section, 2)}_dots_one-sided.png"
    plot_and_save_delaunay(points, tri, title, save_path)

    # Calculate total area
    total_area = calculate_total_area(tri, points)
    print(f"Total Area of {selected_region}_SLIDE{np.round(select_section, 2)}_one-sided:, {total_area} mm2.", file=log_file)

    # Estimate density: calculate volume of a 10 micron slice (volume = length × width × height)
    densities = section['cluster'].value_counts() / (total_area * 0.01)
    df_densities = pd.DataFrame({'cluster': densities.index, 'density': densities.values})

    return df_densities
  
    
def process_slide_side(section, half, side, select_section, selected_region, download_base):
               
    # Filter DataFrame based on x_reconstructed values
    subsection = section[section['x_reconstructed'] <= half] if side == "right_side" else section[section['x_reconstructed'] > half]
    print(f"Number of cells on the {side} of slide {select_section} of {selected_region} region is {subsection.shape[0]}.", file=log_file)
    
    if section.shape[0] > 4: #requires a minimum of 5 points to construct the initial simplex        
        # Convert DataFrame to NumPy array
        points = subsection[['x_reconstructed', 'y_reconstructed', 'z_reconstructed']].values

        # Perform Delaunay triangulation and save image
        perturbed_points = delaunay_triangulation_with_perturbation(points)
        tri = Delaunay(perturbed_points)
        title = f"{selected_region}_SLIDE{np.round(select_section, 2)}_{side}"
        #cluster_as_folder = selected_region.translate(str.maketrans({"/":  r"", "-":  r"", " ":  r"", ",": r""}))
        save_path = f"{download_base}results/density_calculations_1/img/{cluster_as_folder}_{np.round(select_section, 2)}_dots_{side}.png"
        plot_and_save_delaunay(points, tri, title, save_path)

        # Calculate total area
        total_area = calculate_total_area(tri, points)
        print(f"Total Area of {selected_region}_SLIDE{np.round(select_section, 2)}_{side}:, {total_area} mm2.", file=log_file)

        # Estimate density: calculate volume of a 10 micron slice (volume = length × width × height)
        densities = subsection['cluster'].value_counts() / (total_area * 0.01)
        df_densities = pd.DataFrame({'cluster': densities.index, 'density': densities.values})
    else:
        print(f"{selected_region}, {select_section}: {side} side and has not enough cells {subsection.shape[0]}, densities on this side will be 0.", file=log_file)
        densities = section['cluster'].value_counts() * 0 #Since the area is too small, densities will be 0.
        df_densities = pd.DataFrame({'cluster': densities.index, 'density': densities.values})
    
    return df_densities

def process_slide_side_logless(section, half, side, select_section, selected_region, download_base):
               
    # Filter DataFrame based on x_reconstructed values
    subsection = section[section['x_reconstructed'] <= half] if side == "right_side" else section[section['x_reconstructed'] > half]
    print(f"Number of cells on the {side} of slide {select_section} of {selected_region} region is {subsection.shape[0]}.")
    
    if section.shape[0] > 4: #requires a minimum of 5 points to construct the initial simplex        
        # Convert DataFrame to NumPy array
        points = subsection[['x_reconstructed', 'y_reconstructed', 'z_reconstructed']].values

        # Perform Delaunay triangulation and save image
        perturbed_points = delaunay_triangulation_with_perturbation(points)
        tri = Delaunay(perturbed_points)
        title = f"{selected_region}_SLIDE{np.round(select_section, 2)}_{side}"
        #cluster_as_folder = selected_region.translate(str.maketrans({"/":  r"", "-":  r"", " ":  r"", ",": r""}))
        save_path = f"{download_base}results/density_calculations_1/img/{cluster_as_folder}_{np.round(select_section, 2)}_dots_{side}.png"
        plot_and_save_delaunay(points, tri, title, save_path)

        # Calculate total area
        total_area = calculate_total_area(tri, points)
        print(f"Total Area of {selected_region}_SLIDE{np.round(select_section, 2)}_{side}:, {total_area} mm2.")

        # Estimate density: calculate volume of a 10 micron slice (volume = length × width × height)
        densities = subsection['cluster'].value_counts() / (total_area * 0.01)
        df_densities = pd.DataFrame({'cluster': densities.index, 'density': densities.values})
    else:
        print(f"{selected_region}, {select_section}: {side} side and has not enough cells {subsection.shape[0]}, densities on this side will be 0.")
        densities = section['cluster'].value_counts() * 0 #Since the area is too small, densities will be 0.
        df_densities = pd.DataFrame({'cluster': densities.index, 'density': densities.values})
    
    return df_densities


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


# def read_and_test_csv_files(filenames, unique_prefixes, folder_path):
#     result_dataframes = {}
    
#     for prefix in unique_prefixes:
#         matching_files = [filename for filename in filenames if filename.startswith(prefix)]
#         if matching_files:
#             for filename in matching_files:
#                 df = pd.read_csv(os.path.join(folder_path, filename), index_col='cluster')
#                 if df.empty:
#                     print(f"The DataFrame from {filename} is empty")
                
def read_and_test_csv_files(filenames, unique_prefixes, folder_path):
    '''Test all files created above, whether they are empty and the columns are in order.'''
    for prefix in unique_prefixes:
        matching_files = [filename for filename in filenames if filename.startswith(prefix + '_')]
        if matching_files:
            for filename in matching_files:
                df = pd.read_csv(os.path.join(folder_path, filename), index_col='cluster')
                if df.empty:
                    print(f"The DataFrame from {filename} is empty")
                    continue  # Skip further checks if the DataFrame is empty
                
                # Check if 'cluster' is the index and 'density_mm3' is a column
                if df.index.name != 'cluster' or 'density_mm3' not in df.columns:
                    print(f"The DataFrame from {filename} has incorrect index or missing 'density_mm3' column")
                    continue  # Skip further checks if the index or column is incorrect
               
                # Check data types of 'cluster' and 'density_mm3' columns
                if df.index.dtype != object or not pd.api.types.is_numeric_dtype(df['density_mm3']):
                    print(f"The DataFrame from {filename} has incorrect data types")
                    continue  # Skip further checks if data types are incorrect
                
                # If all checks pass, store the DataFrame
                #result_dataframes[filename] = df
                
    print("If there are no print messages before this, all DataFrames / csv files contains cells and densities")

def read_csv_files_into_single_dataframe(filenames, unique_prefixes, folder_path):
    result_dataframes = {}
    
    for prefix in unique_prefixes:
        matching_files = [filename for filename in filenames if filename.startswith(prefix + '_')]
        if matching_files:
            dfs = [pd.read_csv(os.path.join(folder_path, filename), index_col='cluster') for filename in matching_files]
            result_dataframes[prefix] = pd.concat(dfs, ignore_index=True)
    
    return result_dataframes

