#!/bin/bash
#SBATCH --job-name=I0Complete    # Job name      
#SBATCH --nodes=1                        # Run all processes on a single node.  
#SBATCH --array=[0-43]%153      # Run an array job, but with a limitation of 50 arrays x time (to not overload Tigerfish with only one user) \\ LastSendedJob=0 MaxJobsToSend=ceil(82*52/1000)=5
#SBATCH --cpus-per-task=1                # Number of CPU cores per task (multithreaded tasks)
#SBATCH  --output=/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/output/I0slurm-%A_%a.out

python -u CompleteI0.py $SLURM_ARRAY_TASK_ID 0
