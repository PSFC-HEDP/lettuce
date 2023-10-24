#!/bin/bash

module purge
module load anaconda3/2023.07-2
module load slurm/current

PYTHONPATH=. python3 python/start_run.py IRIS "$@"

./check_on_runs.sh
