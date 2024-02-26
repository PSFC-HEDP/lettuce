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
That's fine.

If you're running locally you only need to do the twoth line.
~~~bash
pip install -r requirements.txt
~~~

You'll need the LILAC input deck template.
LILAC input decks aren't export controlled, but LLE still doesn't like me sharing them freely.
If you're using this repository, I assume you have LILAC access;
ask me for the template file and I'll send it to you.
It goes in `lettuce/resources`.

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

## some notes on LILAC and IRIS (or: what the LILAC user guide doesn't want you to know)

Do not call LILAC with a command-line argument!  It will break all of the filenames.

Do not set both `gas_press` and `mat_dens` in a single layer!  LILAC will not respect either of them.

Do not use `zoning='user'` with `feather=.false.`!  This will cause bizarre spline zoning that results in early termination.
Setting `feather=.true.`, even if `xfactor=1`, results in correct equal-mass zoning (or equal-thickness zoning if `eqrad=.true.`).
Note that only the shell needs to have `feather=.true.`.  For some reason it's fine if `feather=.false.` in the fill.
No, I don't know what "feathering" is so please stop asking.

Not all materials are listed in the user guide.
You can get an idea for what materials exist by looking at the opacity tables in `/lle/data/opacity_tables`
(*not* `/theory_codes/opacity_tables` as the user guide suggests),
but note that those won't give you accurate mass numbers or default densities.
The complete list of D³He materials is given, to the best of my ability, in [`material.py`](python/material.py).

The user guide gives a few ways to make custom materials.
The fancy way is to simply pass multiple material codes in a `&prof` namelist and specify their abundances with `mixfrac`.
Unfortunately, this only works if both materials are compatible with `iopac=1` (i.e. have an astrophysical table)
or both are compatible with `iopac=7` (i.e. are a valid input to aplmix).
Also, this mixing routine requires aplmix regardless and I don't think I have access to aplmix.
The more reliable way is to use a negative material code and the `&mater` namelist.
If in a `&prof` namelist the material code is negative (`-1` conventionally) *or* if `opgrp=user`,
then LILAC will expect to see a `&mater` namelist immediately after it.
This namelist is used to provide additional material information,
tho be careful as you need to provide *all* of that information.
Specificly, you need to provide the material nuclide species and abundances,
a positive material code (unclear how this is used),
the full opacity table filename (even if this is inferrable from `iopac` and the material code),
and the list of radiation energy bins for that opacity table
(even tho the opacity table specifies the energy bins and it's the same 48 bins that all of the opacity tables use
(the radiation energy bins in each layer must be the same or LILAC will throw a (surprisingly helpful) error message)).

It wasn't working.  it's doing the correct amount of T but it's filling the rest with D.  I switched the order to put Ti last but no dies.  I added a matname and removed matcod from mater.  maybe next I'll try removing ptrit.  or messing with the matcod in &prof.  ope, setting matcod to -1 made it pure D.  wait, am I not specifying a matcod in &mater?  I put back the matcod in &mater and set matcod in &prof to -1, but agen it's oops all D.  which is weird because I was getting T at least before.
just to be sure, let me now remove matname.
that should set me back to one iteration ago.
tf why is it still hydrogen?  is it the matcod=-1?  is that what screws it up?
maybe I need to stop using the FPOT table and use a different one instead?

The IRIS documentation doesn't state what units the inputs must be in.
In reality they're all SI (velocity in m/s, density in kg/m³, and cetera).
Note that temperature is in J, not K or keV.
