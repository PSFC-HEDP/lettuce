#!/bin/bash

module purge
module load slurm/current

export SLURM_TIME_FORMAT=relative
echo "Code  Name     Slurm_ID     State      Reason     Time     Start       End"
echo "-------------- -------- --------- ----------- -------- --------- ---------"
squeue --user="$(whoami)" --states="all" --sort="i" --format="%14j %.8A %.9T %.11r %.8M %.9S %.9e" | grep "LILAC"
