#!/bin/bash
#SBATCH --job-name=simVarDrive     # Job name      
#SBATCH --nodes=1                        # Run all processes on a single node.  
#SBATCH --array=[0-1000]%50   # Run an array job, but with a limitation of 50 arrays x time (to not overload Tigerfish with only one user) \\ LastSendedJob= JobsToSend=ceil(12*10*4*62/1000)=30
#SBATCH --cpus-per-task=1                # Number of CPU cores per task (multithreaded tasks)
#SBATCH  --output=/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/output2/Simuslurm-%A_%a.out

srun python3 -u MainVarDrive.py $SLURM_ARRAY_TASK_ID 0
