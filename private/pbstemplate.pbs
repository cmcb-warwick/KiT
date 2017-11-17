#!/bin/bash
#PBS -l walltime=2:00:00
#PBS -l pmem=2gb
#PBS -l ncpus=1
#PBS -V
#PBS -e localhost:$HOME/$PBS_JOBNAME.err
#PBS -o localhost:$HOME/$PBS_JOBNAME.out

# Random delay to avoid overloading resources.
sleep $[ ( $RANDOM % 10 )  + 1 ]s

# PBS template for KiT. Edit to suit your configuration.

cd $HOME/kit
matlab -nodesktop -nosplash -singleCompThread -r "js=kitLoadJobset('$JOBSET_FILE');kitRunJobs(js,'subset',$PBS_ARRAYID);"

