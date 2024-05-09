#!/bin/bash

module purge
module load lotus/2.6.3
module load slurm/current

python3 python/start_run.py IRIS "$@"
exit_code=$?
if [ ${exit_code} -ne 0 ]; then
	exit ${exit_code}
fi

./check_on_runs.sh
