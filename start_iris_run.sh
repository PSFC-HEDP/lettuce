#!/bin/bash

module purge
module load lotus/2.6.3
module load slurm/current

python3 python/start_run.py IRIS "$@"

./check_on_runs.sh
