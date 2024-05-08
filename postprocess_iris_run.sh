#!/bin/bash

module purge
module load anaconda3/2020.11b

python3 python/postprocess_iris_run.py "$@"
