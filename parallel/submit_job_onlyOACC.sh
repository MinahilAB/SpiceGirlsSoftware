#!/bin/bash
#SBATCH --job-name=ATM_Model_onlyOACC
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg
#SBATCH --gres=gpu:1
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=output.out
#SBATCH --error=error.err
#SBATCH --account=ICT25_MHPC_0

set -euo pipefail

# Add the path to the SpiceGirlsSoftware directory
SPG_DIR=${HOME}/SpiceGirlsSoftware/parallel

# Where ou want the output to go (can leave unchanged)
OUTDIR=${HOME}/Jobs/ATM_Model_onlyOACC
mkdir -p ${OUTDIR}

# Set nx gtidsize and the simulation time
NX_SIZE=100
SIM_TIME=1000.0

# Set Makefile flags
DEBUG=1
USE_OPENACC=1
USE_OPENMP=0

# Set the number of OMP threads for CPU thread parallelisation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# ===== OpenACC DEBUGGING FLAGS FOR GCC =====
export GOMP_DEBUG=1              # Enable libgomp debugging
export GCC_OFFLOAD_DEBUG=1       # Show offload operations
export CUDA_LAUNCH_BLOCKING=1    # Synchronous kernel launches
# ====================================

echo " ==== Loading modules... ===== "
module purge
module load cuda/12.6
module load gcc/12.2.0
module load netcdf-fortran/4.6.1--openmpi--4.1.6--gcc--12.2.0-spack0.22

echo " ==== Building Project... ===== "
WORK_ROOT=$(mktemp -d "${SCRATCH}/${SLURM_JOB_NAME}_XXXXXX")
mkdir -p "${WORK_ROOT}"
cp -rf "${SPG_DIR}" "${WORK_ROOT}/"

WORK_DIR="${WORK_ROOT}/parallel"
make -C "${WORK_DIR}" clean
make -C "${WORK_DIR}" DEBUG=${DEBUG} USE_OPENACC=${USE_OPENACC} USE_OPENMP=${USE_OPENMP}

srun nsys profile \
    --trace=cuda,nvtx \
    -o "${OUTDIR}/%q{SLURM_JOB_ID}_N%q{SLURM_JOB_NUM_NODES}_%q{SLURM_PROCID}" ${WORK_DIR}/model ${NX_SIZE} ${SIM_TIME}

mv output.nc ${OUTDIR}/${SLURM_JOB_ID}_output.nc
mv statistics_*.txt ${OUTDIR}/${SLURM_JOB_ID}_statistics_*.txt

rm -rf "${WORK_ROOT}"
echo "==== Cleanup Done! ==== "

##SBATCH --gres=gpu:1
##SBATCH --partition=boost_usr_prod
##SBATCH --qos=boost_qos_dbg
##SBATCH --account=ICT25_MHPC_0
#OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}