#!/bin/bash

module purge
module load anaconda3/2023.07-2
module load slurm/current

python3 python/start_run.py LILAC "$@"

./check_on_runs.sh
