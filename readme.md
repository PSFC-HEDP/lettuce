# Lettuce

This is a suite of Bash and Python scripts for running LILAC and IRIS from the command line.

The name is a play on Varchas Gopalaswamy's Lotus
(which is itself a continuation of the theme of LILAC and IRIS),
which has some similar functionality.
Unlike Lettuce, tho, Lotus is a Python *package*.
So it helps you write your own Python scripts, but you still have to write the scripts yourself.
Props to Varchas, because it's a nice and feature-rich package,
but I'm just not about that life, you know?

## File structure

- `data/`  
  All of the simulation inputs and outputs, organized into subfolders by shot and code
- `template/`  
  Templates used to write the batch scripts for each simulation
- `python/`
  Python scripts that get called by the bash scripts in the root directory
- `run_inputs.csv`  
  Table in which the user specifies the parameters for every simulation
- `run_outputs.csv`  
  Alphabetical list of all simulations, their statuses, and their yields, bang-times, etc
- `runs.log`  
  Chronological list of all simulation starts, completions, and cancellations
- `start_lilac_run.sh`  
  Script to arrange the inputs for a LILAC run in a machine-readable format and queue it up
- `start_all_lilac_runs.sh`  
  Script to run all LILAC simulations that are specified in run_inputs.csv but have not yet been run
- `check_on_lilac_runs.sh`  
  Script to report on the status of all current and recent LILAC runs
- `start_iris_run.sh`  
  Script to arrange the inputs for an IRIS run in a machine-readable format and queue it up
- `check_on_iris_runs.sh`  
  Script to report on the status of all current and recent IRIS runs
- `tune_lilac_parameters.sh`  
  Script to queue up and plot an array of LILAC runs to seek out a specific set of observables

## Installation

## Recommended workflow

It will warn you if the simulation appears to be exactly the same as one that has already been run.
