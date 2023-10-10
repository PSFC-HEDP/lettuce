# Lettuce

[![continuous integration](https://github.com/PSFC-HEDP/lettuce/actions/workflows/main.yml/badge.svg)](https://github.com/PSFC-HEDP/lettuce/actions/workflows/main.yml)

This is a suite of Bash and Python scripts for running LILAC and IRIS from the command line.

The name is a play on Varchas Gopalaswamy's Lotus
(which is itself a continuation of the theme of LILAC and IRIS),
which has some similar functionality.
Unlike Lettuce, tho, Lotus is a Python *package*.
So it helps you write your own Python scripts, but you still have to write the scripts yourself.
Props to Varchas, because it's a nice and feature-rich package,
but I'm just not about that life, you know?

## File structure

- `runs/`  
  All of the simulation inputs and outputs, organized into subfolders by shot and code
- `templates/`  
  Templates used to write the batch scripts for each simulation
- `pulse_shapes/`  
  Pulse shape files downloaded from OmegaOps
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

The first 

It will warn you if the simulation appears to be exactly the same as one that has already been run.

## Run inputs

Runs are defined in run_inputs.csv.
Each row contains a unique name for the run followed by all the information needed to run the simulation.
The columns are as follows:

- **laser energy**  
   The total time-integrated energy incident on the capsule, in kilojoules.
- **pulse shape**  
   The name of the laser pulse shape.  The 1 ns square pulse is "SG10v001".  If you're not sure what the name of your pulse shape is, consult the [OMEGA pulse shape library](https://omegaops.lle.rochester.edu/cgi-script/pulseShapes).
- **outer diameter**  
  The size of the capsule, in micrometers.
- **shell material**  
  The name of the shell material.  Many materials have multiple acceptable names (for example, "SiO2" and "glass" both do the same thing).
- **shell thickness**  
  The thickness of the shell, in micrometers.
- **aluminum thickness**  
  The thickness of the aluminum coating on the outside of the shell, in micrometers (usually 0.1).
- **fill**  
  A string that specifies the density and composition of the gas fill.  It should look something like this: "12atm 3He + 6atm D".  Each component is given as a molecular pressure at 293K followed by the name of the element.  Note that these are molecular pressures.  "D" and "D2" are interchangeable.
- **absorption fraction**  
  The ratio between the laser energy absorbed by the capsule and the total laser energy.
- **flux limiter**  
  The free-streaming electron sharp-cutoff flux limiter coefficient.  Units unknown.  Setting it to 0 will use Valeri Goncharov's nonlocal model instead.
- **laser degradation**  
  If this is nonzero, the laser pulse will be cut short by this many picoseconds.  This is primarily used for matching simulated yields to experimental ones.
