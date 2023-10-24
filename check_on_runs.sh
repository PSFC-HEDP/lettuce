#!/bin/bash

module purge
module load slurm/current

export SLURM_TIME_FORMAT=relative
queue="$(squeue --user="$(whoami)" --states=all --sort=i --format='%14j %.8A %.9T %.11r %.8M %.9S %.9e')"
echo "Code  Name     Slurm_ID     State      Reason     Time     Start       End"
echo "-------------- -------- --------- ----------- -------- --------- ---------"
grep LILAC <<< "$queue"
grep IRIS <<< "$queue"
