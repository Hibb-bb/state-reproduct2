#!/bin/bash
#SBATCH --account=p32234
#SBATCH --job-name=cpa
#SBATCH --nodes=1
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=slurm-cpa.out
#SBATCH --error=slurm-cpa.err

module purge all
module load gcc/11.2.0
source ~/.bashrc

cd /projects/b1094/ywl7940/state-reproduce
cd baselines

source ./.venv/bin/activate


HYDRA_FULL_ERROR=1 sh scripts/train.sh cpa replogle 1