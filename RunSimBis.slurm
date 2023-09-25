#!/bin/bash
#SBATCH --job-name=sim2New     # Job name      
#SBATCH --nodes=1                        # Run all processes on a single node.  
#SBATCH --array=[0-50]%50    # Run an array job, but with a limitation of 50 arrays x time (to not overload Tigerfish with only one user) \\ LastSendedJob=3 JobsToSend=ceil(82*52*6/1000)=26
#SBATCH --cpus-per-task=1                # Number of CPU cores per task (multithreaded tasks)
#SBATCH  --output=/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/output/Simuslurm-%A_%a.out

srun python -u MainBis.py $SLURM_ARRAY_TASK_ID 0
