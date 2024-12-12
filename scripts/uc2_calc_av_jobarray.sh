#!/bin/bash -l

#SBATCH --job-name=mtrx_crunch
#SBATCH --account=proj134
#SBATCH --partition=prod
#SBATCH --constraint=uc2
#SBATCH --mem=0 
#SBATCH --time=24:00:00
#SBATCH --output=outputs/%A_%a.log
#SBATCH --error=error/stderr-%A_%a.log
#SBATCH --array=0-4%5

echo "Copying files now.."

###
# cp files to be read by the code to SHMDIR to avoid overloading the gpfs by read write requests
#provide this SHMDIR to the code with os.environ.get("SHMDIR", "")
#Alternatively, #SBATCH --constraint=nvme will give me $TMPDIR which is slower but has more space
#Each allocated node has its own SHMDIR and/or TMPDIR. only gpfs is shared among allocated nodes.
cd /gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain
cp --parents ./releases/20230830/manifest.json $SHMDIR
cp --parents ./metadata/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv $SHMDIR
cp --parents ./metadata/WMB-10X/20230830/gene.csv $SHMDIR
cp --parents ./expression_matrices/WMB-10Xv2/20230630/*log2* $SHMDIR
cp --parents ./expression_matrices/WMB-10Xv3/20230630/*log2* $SHMDIR
cp --parents ./expression_matrices/WMB-10XMulti/20230830/*log2* $SHMDIR
mkdir -p $SHMDIR/metadata/averages/mean
mkdir -p $SHMDIR/metadata/averages/median
###

echo "Finished copying files, loading Python.."

# Load Python environment
module unload python #because of my default .zshrc 
conda activate cell2loc

#Provide cycles (where to parse the data) which is passed on the the working node
cycles=( 1163 1165 1167 1184 1186)
celltype_batch=${cycles[${SLURM_ARRAY_TASK_ID}]}
echo $celltype_batch

# python /gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/notebooks/scripts/calc_av_jobarray_writegpfs.py "$celltype_batch"

python /gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/notebooks/scripts/calc_av_jobarray_classesonly_writegpfs.py "$celltype_batch"

#After the code has run on the nodes, CP results to gpfs 
cp -R $SHMDIR/metadata/averages /gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/metadata/averages/20230830/

# Print "Job's done" message to the log files
echo "All files were copied to gpfs!"