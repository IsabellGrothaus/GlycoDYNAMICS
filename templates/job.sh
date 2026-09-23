#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --account=hbb00004

unset SLURM_EXPORT_ENV

module load gcc/11.3.0 cuda/11.8 openmpi/gcc.11/4.1.4
module load gromacs/2022.5

export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps.$SLURM_JOB_ID
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log.$SLURM_JOB_ID

nvidia-cuda-mps-control -d

export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_GPU_FORCE_UPDATE_DEFAULT_GPU=true

export OMP_NUM_THREADS=6

mpiexec -np 12 gmx_mpi_plumed mdrun -ntomp $OMP_NUM_THREADS -s topol.tpr -pme gpu -nb gpu -update cpu -bonded gpu -pin on -pinstride 1 -plumed plumed_allangle.dat -multidir TOPO? TOPO?? -replex 400 -hrex -dlb no -maxh 24.0

echo quit | nvidia-cuda-mps-control


