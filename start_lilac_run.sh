#!/bin/bash

module purge
module load anaconda3/2020.11b
module load slurm/current

python3 python/start_run.py LILAC "$@"
exit_code=$?
if [ ${exit_code} -ne 0 ]; then
	exit ${exit_code}
fi

./check_on_runs.sh
