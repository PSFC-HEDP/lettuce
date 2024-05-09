# Lettuce

[![test status](https://github.com/PSFC-HEDP/lettuce/actions/workflows/main.yml/badge.svg)](https://github.com/PSFC-HEDP/lettuce/actions/workflows/main.yml)

This is a suite of Bash and Python scripts for running LILAC and IRIS from the command line.

The name is a play on Varchas Gopalaswamy's Lotus
(which is itself a continuation of the theme of LILAC and IRIS),
which has some similar functionality.
Unlike Lettuce, tho, Lotus is a Python *package*.
So it helps you write your own Python scripts, but you still have to write the scripts yourself.
Not so with Lettuce!

## Repository contents

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
Right now I'm working on `anaconda3/2020.11b`, so that would look like this:
~~~bash
module load anaconda3/2020.11b
python3 -m pip install -r requirements.txt
~~~
It will probably say something about user installation because normal site-packages is not writeable.
That's fine.

Occasionally LLE will have installed incompatible versions of two packages on a single module.
When that happens, you can either take it up with Jonathan Carrol-Nellenback or start throwing around some `--upgrade`s and hope it helps.
In my case, I had to `pip install --upgrade h5py` to make it play nice with NumPy.

If you're running locally you only need to do the twoth line.
~~~bash
pip install -r requirements.txt
~~~

You'll need the LILAC input deck template.
LILAC input decks aren't export controlled, but LLE still doesn't like me sharing them freely.
If you're using this repository, I assume you have LILAC access;
ask me for the template file and I'll send it to you.
It goes in `lettuce/resources/templates/`.

## Recommended workflow

### LILAC

The first step is to describe the run you want to do in `run_inputs.csv`.
See [§ Run inputs](#run-inputs) for more information.

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
PYTHONPATH=python/ python -m pytest
~~~
On Windows PowerShell the command is
~~~powershell
$env:pythonpath="python\"; python -m pytest
~~~

### Directly editing input decks

In the case of both LILAC and IRIS, you should ideally be able to specify all of the options you need from the command line or using `run_inputs.csv`.
However, LILAC and IRIS are far more featureful than my code, and thus you may occasionally need to edit the input decks directly.
All you'll need to do is go into the relevant subdirectory of `runs/`,
find `lilac_input_deck.txt` in the case of LILAC or `inputdeck.txt` in the case of IRIS,
make the desired change,
and then call `sbatch run.sh`.
That script contains all of the slurm options, the updates to `runs.log` and `run_outputs.csv`, and the postprocessing step.

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
  The state of the most recent Slurm job, followed by either "LILAC" or "IRIS" depending on what that job was.
  States are one of "pending", "running", "completed", "failed", or "timed out".
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

## To do

I have zero time to work on this rite now, but there are a few changes and enhancements I want to make that will greatly improve this library before I hand LILAC off.  and some major changes.
- allow multiple runs be started simultaneusly by simply separating their names with commas.
- have options like `--flux_limiter=.07` to automaticly generate a new set of run_inputs based on an existing row in run_inputs.csv but with one parameter tweaked and an appropriate tag appended to the name with a slash.
- make the output PDFs portait instead of landscape, or use two columns.
- write an IRIS postprocessing script that makes a PDF with some spectra and images.
- change the status recorded in run_outputs.csv from "pending" to "running" when it gets off the queue (right now it goes straight from "pending" to "completed" when it finishes, which is misleading.)
- find a way to change the status recorded in run_outputs.csv to "cancelled" when it gets cancelled (it's rough because none of my code gets called when that happens).

## some notes on LILAC and IRIS (or: one weird trick the LILAC user guide doesn't want you to know)

Do not call LILAC with a command-line argument!  It will break all of the filenames.

Do not set both `gas_press` and `mat_dens` in a single layer!  LILAC will not respect either of them.

Do not use `zoning='user'` with `feather=.false.`!  This will cause bizarre spline zoning that results in early termination.
Setting `feather=.true.`, even if `xfactor=1`, results in correct equal-mass zoning (or equal-thickness zoning if `eqrad=.true.`).
Note that only the shell needs to have `feather=.true.`.  For some reason it's fine if `feather=.false.` in the fill.
No, I don't know what "feathering" is so please stop asking.

Sometimes a run will terminate before shock breakout, so you'll see the outer layer ablating for 800 ps or so but noting else.
This is especially common when the shell has multiple layers,
and especially especially when at least one of the layers is strong CD.
I know it has to do with the density and number of zones in each layer,
but I haven't found a way to reliably predict or prevent it.

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
Note that temperature must be in J, not K or keV.

IRIS generates a HDF5 file with two main datasets that give all of the imaging and spectral information.
The dataset `/images/image` contains all of the images.
Each value is the space- and energy-resolved intensity at a given point at a given time measured by a given detector, measured in particles/MeV/m^2/sr.
The values are arranged in 6 dimensions as follows:
- dimension 0 indexes the different times of the implosion (given in nanoseconds in `/time`);
- dimension 1 indexes the different reaction products (0: primary fusion neutrons, 1: proton-scattered neutrons, 2: deuteron-scattered neutrons, 3: triton-scattered neutrons, 4: carbon-scatterd neutrons, 5: neutrons from deuteron breakup, 6: knock-on deuterons, 7: knock-on tritons);
- dimension 2 indexes the different detection energy bins (the bin edges are given in MeV in `/images/energy`);
- dimension 3 indexes the different detectors (the position of each is given in radians in `/images/theta` and `/images/phi`);
- dimension 4 indexes the different x bins (the edges are evenly spaced and centered on x = 0; the absolute value of the outermost bin edges is given in meters in `/images/size`);
- dimension 5 indexes the different y bins (the edges are evenly spaced and centered on y = 0; the absolute value of the outermost bin edges is given in meters in `/images/size`).

The dataset `/spectra/dNdE` contains all of the spectra.
Each value is the energy-resolved fluence of a given particle species at a given time measured by a given detector, in unknown units.
The values are arranged in 6 dimensions as follows:
- dimension 0 indexes the different times of the implosion (see above);
- dimension 1 indexes the different fusion reactions (0 is DT, 1 is DD, I don't remember what 2 is);
- dimension 2 indexes particles by the number of scattering events;
- dimension 3 indexes the different types of reaction product (see above);
- dimension 4 indexes the different energy bins (the bin edges are given in joules in `/spectra/energy`);
- dimension 5 indexes the different detectors (the position and size of each is given in `/spectra/theta`, `/spectra/phi`, `/spectra/min_solid_angle`, and `/spectra/max_solid_angle`).

In addition to spectra, IRIS also returns the average ρL seen by each detector as the dataset `/spectra/rhoL`.
Each value is the average line-integrated density experienced by particles collected by a given detector, at a given time.
The values are arranged in 3 dimensions as follows:
- dimension 0 indexes the different times of the implosion (see above);
- dimension 1, based on its length, must either index by fusion reaction or by number of scattering events, but since the output file isn't annotated it's impossible to tell;
- dimension 2 indexes the different detectors (see above).
