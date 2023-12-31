# Lettuce

[![test status](https://github.com/PSFC-HEDP/lettuce/actions/workflows/main.yml/badge.svg)](https://github.com/PSFC-HEDP/lettuce/actions/workflows/main.yml)

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

The only dependencies are on PyPI.
If you're running on Bluehive, you should make sure you have all of the dependencies installed on the relevant module.
Right now I'm working on `anaconda3/2023.07-2`, so that would look like this:
~~~bash
module load anaconda3/2023.07-2
python3 -m pip install -r requirements.txt
~~~~
It will probably say something about user installation because normal site-packages is not writeable.

If you're running locally you only need to do the twoth line.
~~~bash
pip install -r requirements.txt
~~~

## Recommended workflow

### LILAC

The first step is to describe the run you want to do in `run_inputs.csv`.
See [§ Run inputs](#Run inputs) for more information.

To start the run, call
~~~bash
./start_lilac_run.sh [NAME]
~~~
The `[NAME]` should be the same as what's in the "name" column of `run_inputs.csv`.
It can also end in an asterisk to run all pattern-matched names in the table.
This currently only works for a single asterisk at the end of the name.
It will warn you if the simulation appears to be exactly the same as one that has already been run.
If you have multiple runs you want to do, call this once for each one and they will be queued in parallel.

While you're waiting, you can see how your runs are doing with
~~~bash
./check_on_runs.sh
~~~
This will print out all of the LILAC that is currently queued or running.

Once the LILAC job finishes, it will automatically post-process the result
and generate a PDF containing numbers of interest in its run directory.

### IRIS

After LILAC is done, you can run IRIS.
Eventually it would be nice if it automatically ran IRIS after LILAC, but for now it's manual.
To start an IRIS run, call
~~~bash
./start_iris_run.sh [NAME]
~~~

While you're waiting, you can see how your runs are doing with
~~~bash
./check_on_runs.sh
~~~
This will print out all of the IRIS that is currently queued or running.

Once the IRIS job finishes, it will automatically post-process the result
and generate a PDF containing spectra and images in its run directory.

### Unit tests

If you want to run unit tests, mind the PythonPath.  On Unix the command is
~~~bash
PYTHONPATH=python python -m pytest
~~~
On Windows PowerShell the command is
~~~powershell
$env:pythonpath="python"; python -m pytest
~~~

## Run inputs

Runs are defined in `run_inputs.csv`.
Each row contains a unique name for the run followed by all the information needed to run the simulation.
The columns are as follows:

- **name**  
  A uniquely identifying string for this simulation.
- **laser energy**  
  The total time-integrated energy incident on the capsule, in kilojoules.
- **pulse shape**  
  The name of the laser pulse shape.  The 1 ns square pulse is "SG10v001".  If you're not sure what the name of your pulse shape is, consult the [OMEGA pulse shape library](https://omegaops.lle.rochester.edu/cgi-script/pulseShapes).
- **beam profile**  
  The name of the beam profile.  For most people it will be "SG5 SSD".
- **outer diameter**  
  The size of the capsule, in micrometers.
- **shell material**  
  The name or code of the shell material.  Some materials have multiple acceptable names (for example, "SiO2" and "glass" both do the same thing).
  For layered shells, you can put multiple materials separated by slashes (for example, "CH/148/CH").
- **shell thickness**  
  The thickness of the shell, in micrometers.  If multiple materials were given, the same number of thicknesses must be given, separated by slashes.
- **aluminum thickness**  
  The thickness of the aluminum coating on the outside of the shell, in micrometers (usually 0.1).
- **fill**  
  A string that specifies the density and composition of the gas fill.  It should look something like this: "12atm 3He + 6atm D".  Each component is given as a molecular pressure at 293K followed by the name of the element.  Note that these are molecular pressures.  "D" and "D2" are interchangeable.
- **absorption fraction**  
  The ratio between the laser energy absorbed by the capsule and the total laser energy.
- **flux limiter**  
  The optional free-streaming electron sharp-cutoff flux limiter coefficient.  Units unknown.  Setting it to 0 will use Valeri Goncharov's nonlocal model instead.
- **laser degradation**  
  If this is given and nonzero, the laser pulse will be cut short by this factor.  Inputting 1 will completely turn off the laser.  This is primarily used for matching simulated yields to experimental ones.
- **shell density multiplier**  
  The optional factor by which the shell density should differ from whatever it normally is for that material.  This can be used to match simulated ρRs to experimental ones.

## Key run outputs

The run outputs are described in detail in the `output.pdf` file of each run directory.
For a conglomerated summary, tho, one can look at `run_outputs.csv`.
The columns are as follows:
- **name**  
  The uniquely identifying string for this simulation.
- **slurm ID**  
  The ID number assigned by Slurm to the most recent job that was submitted related to this.
- **status**  
  The state of the Slurm job; one of "pending", "running", "completed", "failed", or "timed out".
- **status changed**  
  The date and time at which the status of the Slurm job last changed.
- **yield**  
  The nuclear yield of the implosion (for whichever nuclear reaction had the highest yield).
- **bang-time**  
  The time of peak nuclear emission (for whichever nuclear reaction had the highest yield), in nanoseconds.
- **convergence ratio**  
  The ratio between the initial radius of the gas fill and the minimum radius of the gas fill.
- **areal density**  
  The burn-averaged ρR, in milligrams per square centimeter.
- **ion temperature**  
  The burn-averaged ion temperature, in kiloelectron-volts.

## some notes on LILAC and IRIS

Do not call LILAC with a command-line argument!  It will break all of the filenames.

Do not set both `gas_press` and `mat_dens` in a single layer!  LILAC will not respect either of them.

Do not use `zoning='user'` with `feather=.false.`!  This will cause bizarre spline zoning that results in early termination.
Setting `feather=.true.`, even if `xfactor=1`, results in correct equal-mass zoning (or equal-thickness zoning if `eqrad=.true.`).
Note that only the shell needs to have `feather=.true.`.  For some reason it's fine if `feather=.false.` in the fill.
No, I don't know what "feathering" is so please stop asking.

Not all materials are listed in the user guide.
You can get an idea for what materials exist by looking at the opacity tables in `/lle/data/opacity_tables`,
but note that those won't give you accurate mass numbers or default densities.
The complete list of D³He materials is given, to the best of my ability, in [`material.py`](python/material.py).

The user guide gives a few ways to make custom materials, none of which work.
The real way to make a custom material is

The IRIS documentation doesn't state what units the inputs must be in.
In reality they're all SI (velocity in m/s, density in kg/m³, and cetera).
Note that temperature is in J, not K.
