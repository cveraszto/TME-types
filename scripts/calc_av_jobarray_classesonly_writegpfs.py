import sys
import json
import os
import pandas as pd
import anndata
import time
import numpy as np
from scipy import stats

download_base = os.environ.get("SHMDIR", "/gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/")

write = "/gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/"

url = download_base + '/releases/20230830/manifest.json'
with open(url, 'r') as json_file:
    manifest = json.load(json_file)

metadata = manifest['file_listing']['WMB-10X']['metadata']

rpath = metadata['cell_metadata_with_cluster_annotation']['files']['csv']['relative_path']
#This file countains all sc-seq cells and metadata
file = os.path.join( download_base, rpath)
cell = pd.read_csv(file,dtype={"neurotransmitter":str})
cell.set_index('cell_label',inplace=True)
#This is where the raw data is stored
matrices = cell.groupby(['dataset_label','feature_matrix_label'])[['library_label']].count()
matrices.columns  = ['cell_count']
#As a change we will calculate subclass/supertype level averages
#clusters = np.unique(cell['subclass'])
clusters = np.unique(cell['supertype'])
#Full gene set used during sequencing
rpath = metadata['gene']['files']['csv']['relative_path']
file = os.path.join( download_base, rpath)
genes = pd.read_csv(file,dtype={"comment":str})
genes.set_index('gene_identifier',inplace=True)

total_start = time.process_time()
# Check for command-line argument
if len(sys.argv) < 2:
    print("Please provide a value for x.")
else:
    x = int(sys.argv[1]) # Convert x to an integer
    clusters_to_process = clusters[x:x+1]
    
    for cluster in clusters_to_process:

        # Initialize gdata as an empty gene expression DataFrame
        gdata = pd.DataFrame(index=cell[cell['supertype'] == cluster].index, columns=genes.index)
        
        print(f"matrix size for {cluster} is: {gdata.shape}", flush=True) #make sure it is printed as soon as possible (will be in the logs)
        for matindex in matrices.index :

            ds = matindex[0]
            mp = matindex[1]

            print(mp, flush=True)

            expression_matrices = manifest['file_listing'][ds]['expression_matrices']
            rpath = expression_matrices[mp]['log2']['files']['h5ad']['relative_path']
            file = os.path.join( download_base, rpath)

            start = time.process_time()
            ad = anndata.read_h5ad(file,backed='r')
            exp = ad[ :, genes.index ].to_df()
            exp = exp[ exp.index.isin(cell[cell['supertype'] == cluster].index) ]
            gdata.loc[ exp.index, genes.index ] = exp
            print(" - time taken: ", time.process_time() - start, flush=True)

            ad.file.close()
            del ad

        clustername = cluster.translate(str.maketrans({"/":  r"", "-":  r"_", " ":  r"_"}))

        median = gdata.median(axis=0)
        median = median.to_frame(name=cluster)
        file = f"{write}metadata/averages/20230830/supertype_median/{clustername}.csv"
        median.to_csv(file, header=True, index=True)

        trim_means = stats.trim_mean(gdata, 0.25)
        trim_means = pd.DataFrame(trim_means, index=median.index, columns=[cluster])
        file = f"{write}/metadata/averages/20230830/supertype_mean/{clustername}.csv"
        trim_means.to_csv(file, header=True, index=True)

        del median, trim_means, gdata

        print("total time taken: ", time.process_time() - total_start, flush=True)
