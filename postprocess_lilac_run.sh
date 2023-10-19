#!/bin/bash

module purge
module load anaconda3/2023.07-2
module load slurm/current

PYTHONPATH=. python3 python/postprocess_lilac_run.py "$@"
