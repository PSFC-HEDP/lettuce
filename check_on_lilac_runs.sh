#!/bin/bash

module purge
module load slurm/current

squeue --user="$(whoami)" --states="all" --name="LILAC" --Format="JobID:.9,State:.12,Reason:.12,TimeUsed:.9,StartTime:.21,EndTime:.21"
