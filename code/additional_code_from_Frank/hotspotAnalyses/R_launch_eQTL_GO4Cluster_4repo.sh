#!/bin/bash -l
#PBS -l walltime=15:00,nodes=1:ppn=1,mem=1gb
#PBS -o logs/
#PBS -e logs/

# launch this shell for each hotspot:
# for i in {1..102}; do qsub -v j=${i} R_launch_eQTL_GO4Cluster_4repo.sh; done

cd ~/falbert/eQTL1k/scripts

module load R
Rscript eQTL_GO4Cluster_170420_4repo.R $j
wait
