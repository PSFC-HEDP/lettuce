#!/bin/bash

module purge
module load slurm/current

export SLURM_TIME_FORMAT=relative
queue="$(squeue --user="$(whoami)" --states=all --sort=i --format='%.17j %.8A %.9T %.11r %.8M %.10S %.10e')"
lilac_runs="$(grep LILAC <<< "$queue")"
iris_runs="$(grep IRIS <<< "$queue")"
if [ "$lilac_runs" != "" ] || [ "$iris_runs" != "" ]; then
	echo " Code        Name Slurm ID     State      Reason     Time      Start        End"
	echo "-----------------+--------+---------+-----------+--------+----------+----------"
	if [ "$lilac_runs" != "" ]; then
		echo "$lilac_runs"
	fi
	if [ "$iris_runs" != "" ]; then
		echo "$iris_runs"
	fi
else
	echo "No active or recent runs"
fi
