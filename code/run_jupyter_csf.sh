#!/usr/bin/env bash
#SBATCH --job-name="csf_jup"
#SBATCH --chdir=/scratch/devel/pnieto/projects/CSF/
#SBATCH --error=/scratch/devel/pnieto/projects/CSF/output/integration/logs/%x_%J.err
#SBATCH --output=/scratch/devel/pnieto/projects/CSF/output/integration/logs/%x_%J.out
#SBATCH --ntasks=1
#SBATCH -c 48
#SBATCH --mem=250G
#SBATCH --time=12:00:00

source /software/crgadm/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
conda activate /scratch_isilon/groups/singlecell/shared/conda_env/csf

echo [`date "+%Y-%m-%d %T"`] started job on $HOSTNAME

ulimit -n 16000
export HDF5_USE_FILE_LOCKING="FALSE"

echo running script $1

jupyter nbconvert --execute --to html $1 

echo [`date "+%Y-%m-%d %T"`] finished job