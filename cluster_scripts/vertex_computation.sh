#!/bin/bash
#SBATCH -A def-vidotto
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=0-01:00:00
#SBATCH --job-name=BWH_computation
#SBATCH --output=BWH_computation.log
#SBATCH --error=BWH_computation.err

# folders

ROOT_DIR=/home/frisus95/projects/def-vidotto/frisus95
JULIA_DIR=${ROOT_DIR}/julia-1.8.0
SL2CFOAM_DIR=${ROOT_DIR}/sl2cfoam-next-dev
FASTWIG_TABLES_PATH=${SL2CFOAM_DIR}/data_sl2cfoam

export LD_LIBRARY_PATH="${SL2CFOAM_DIR}/lib":$LD_LIBRARY_PATH
export JULIA_LOAD_PATH="${SL2CFOAM_DIR}/julia":$JULIA_LOAD_PATH

setrpaths.sh --path ${JULIA_DIR} [--add_origin]
setrpaths.sh --path ${SL2CFOAM_DIR} [--add_origin]

# otherwise libmpc.so.3 is not found
export LD_LIBRARY_PATH="/cvmfs/soft.computecanada.ca/gentoo/2020/usr/lib":$LD_LIBRARY_PATH

echo "Running on: $SLURM_NODELIST"
echo


# parameters

JULIA_DIR=${JULIA_DIR}
BASE_DIR=${ROOT_DIR}/BlackToWhiteHole
SL2CFOAM_DATA_DIR=${SLURM_TMPDIR}
STORE_FOLDER=${BASE_DIR}
T_PARAMETER=100

CODE_TO_RUN=vertex_computation.jl



# loading fastwig tables

echo "Copying fastwig tables to: $SLURM_TMPDIR ..."
echo

cp ${FASTWIG_TABLES_PATH}/* $SLURM_TMPDIR/



# start computation

echo "Starting computation..."
echo

${JULIA_DIR}/bin/julia -p $SLURM_TASKS_PER_NODE ${BASE_DIR}/bin/${CODE_TO_RUN} $BASE_DIR $SLURM_TMPDIR $T_PARAMETER



# compressing amplitudes

echo "Compressing and copying computed amplitudes to ${BASE_DIR}..."
echo

cd $SLURM_TMPDIR

tar -czvf ${BASE_DIR}/vertex_booster_tensors.tar.gz vertex/



echo "Completed"
echo
