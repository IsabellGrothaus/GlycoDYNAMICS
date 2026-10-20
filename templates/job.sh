#!/bin/bash -l
#SBATCH --gres=gpu:a100:4
#SBATCH --time=24:00:00
#SBATCH --job-name=glycan
#SBATCH --export=NONE

#alex has 16 CPU per GPU
#-np should be number of replicas
#8 gpu per node
#16 cpu per gpu
#1 node per simulation

unset SLURM_EXPORT_ENV

module purge
module load gromacs/2022.5-gcc11.5.0-openmpi4-mkl-cuda-plumed2.8.3

export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps.$SLURM_JOB_ID
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log.$SLURM_JOB_ID

nvidia-cuda-mps-control -d

export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_GPU_FORCE_UPDATE_DEFAULT_GPU=true

#Production
mpirun -np 12 gmx_mpi mdrun -s topol.tpr -pme gpu -nb gpu -update cpu -bonded gpu -pin on -pinstride 1 -plumed plumed_allangles.dat -multidir TOPO? TOPO?? -replex 400 -hrex -dlb no -maxh 24

#Restart
#mpirun -np 12 gmx_mpi mdrun -cpi state.cpt -s topol.tpr -pme gpu -nb gpu -update cpu -bonded gpu -pin on -pinstride 1 -plumed plumed_allangles.dat -multidir TOPO? TOPO?? -replex 400 -hrex -dlb no -maxh 24

# this will stop the MPS-server
echo quit | nvidia-cuda-mps-control


