#!/bin/bash
#SBATCH --job-name=FullFI #Name of the job
#SBATCH --clusters=fisica
#SBATCH --partition=cpu.cecc # Parition (queue) name
#SBATCH --nodes=1 # Number of nodes
#SBATCH --ntasks=1 # Number of tasks (processes)
#SBATCH --cpus-per-task=20 # Number of cpu cores per task
#SBATCH --mem=5G # Memory requiered per node (1GB in this case)
#SBATCH --time=14:01:00 # Time limit (HH::MM:SS)
#SBATCH --output=/homes/fisica/clviviescasr/ninino/%j_%x.out # Standard output file (%j will be replaced by job ID)
#SBATCH --error=/homes/fisica/clviviescasr/ninino/%j_%x.err # Standard error file
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ninino@unal.edu.co
#SBATCH --export=SCRATCH_DIR=/scratch/$SLURM_JOB_ACCOUNT/$SLURM_JOB_USER/$SLURM_JOB_ID 
 cd $SCRATCH_DIR

#first, unload Modules
module purge
# Load required modules
module load lang/julia/1.11.3
export JULIA_DEPOT_PATH="/scratchsan/clviviescasr/ninino/julia"

echo "Job started at: $(date)"
echo "Running on node $(hostname)"
echo "Current Working Directory $(pwd)"

# Your application command
julia Code/sampler.jl 
echo "End time: $(date)"
