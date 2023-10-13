#!/bin/bash
#SBATCH --job-name=IRIS
#SBATCH --time=20:00:00
#SBATCH --partition=blizzard
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4GB
#SBATCH --qos=b-standard
#SBATCH --output=<<directory>>/iris_%A.log
#SBATCH --mail-type=END

# Load IRIS and clean up the environment
module purge
module load hdf5/1.12.0/b3
module load iris/2021.09.03/b1

# move to the directory where runs happen
cd "<<directory>>" || exit 1

# run IRIS
iris
exit_code=$?

# figure out if the run was successful
if [ $exit_code -eq 0 ]; then
	exit_string="IRIS run '<<name>>' completes successfully."
else
	exit_string="IRIS run '<<name>>' fails with error_code ${exit_code}."
fi

# log the result to the slurm log and to runs.log
echo "${exit_string}"
echo "$(date +'%m-%d %H:%M') | ${exit_string}" >> ../../../runs.log

# Move the file out of the "output" directory
if [ $exit_code -eq 0 ]; then
	if [ -e output/iris.hdf5 ]; then
		mv output/iris.hdf5 output.h5
	fi
	rm --recursive --force output
fi

exit ${exit_code}
