#!/bin/bash -l

#SBATCH --job-name=mtrx_crunch
#SBATCH --account=proj134
#SBATCH --partition=prod
#SBATCH --constraint=cpu
#SBATCH --mem=0 
#SBATCH --time=12:00:00
#SBATCH --output=outputs/%A_%a.log
#SBATCH --error=error/stderr-%A_%a.log
#SBATCH --array=0-399%400

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
cycles=( 9 19 29 39 49 59 69 79 89 99 109 119 129 139 149 159 169 179 189 199 209 219 229 239 249 259 269 279 289 299 309 319 329 339 349 359 369 379 389 399 409 419 429 439 449 459 469 479 489 499 509 519 529 539 549 559 569 579 589 599 609 619 629 639 649 659 669 679 689 699 709 719 729 739 749 759 769 779 789 799 809 819 829 839 849 859 869 879 889 899 909 919 929 939 949 959 969 979 989 999 1009 1019 1029 1039 1049 1059 1069 1079 1089 1099 1109 1119 1129 1139 1149 1159 1169 1179 1189 1199  )
celltype_batch=${cycles[${SLURM_ARRAY_TASK_ID}]}
echo $celltype_batch

# python /gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/notebooks/scripts/calc_av_jobarray_writegpfs.py "$celltype_batch"

python /gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/notebooks/scripts/calc_av_jobarray_classesonly_writegpfs.py "$celltype_batch"

#After the code has run on the nodes, CP results to gpfs 
cp -R $SHMDIR/metadata/averages /gpfs/bbp.cscs.ch/data/project/proj84/csaba/aibs_10x_mouse_wholebrain/metadata/averages/20230830/

# Print "Job's done" message to the log files
echo "All files were copied to gpfs!"