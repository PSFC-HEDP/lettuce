#!/bin/bash
#SBATCH --job-name=IRIS_<<name>>
#SBATCH --time=10:00:00
#SBATCH --partition=blizzard
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4GB
#SBATCH --qos=b-standard
#SBATCH --output=<<root>>/runs/<<name>>/<<folder>>/iris_%A.log
#SBATCH --mail-type=END

# Load IRIS and clean up the environment
module purge
module load hdf5/1.12.0/b3
module load iris/2021.09.03/b1

# move to the directory where runs happen
cd "<<root>>/runs/<<name>>/<<folder>>" || exit 1

echo "IRIS run '<<name>>' starts."
echo "$(date +'%m-%d %H:%M') | IRIS run '<<name>>' starts." >> "<<root>>/runs.log"

# run IRIS
timeout 9h iris
exit_code=$?

# figure out if the run was successful
if [ $exit_code -eq 0 ]; then
	result="completed"
	exit_string="IRIS run '<<name>>' completes successfully."
elif [ $exit_code -eq 124 ]; then
	result="timeout"
	exit_string="IRIS run '<<name>>' times out."
elif [ $exit_code -eq 137 ]; then
	result="cancelled"
	exit_string="IRIS run '<<name>>' is cancelled."
else
	result="failed"
	exit_string="IRIS run '<<name>>' fails with error_code ${exit_code} (see '<<root>>/runs/<<name>>/<<folder>>/iris_$SLURM_JOB_ID.log')."
fi

# log the result to the slurm log and to runs.log
echo "${exit_string}"
echo "$(date +'%m-%d %H:%M') | ${exit_string}" >> "<<root>>/runs.log"

# Move the file out of the "output" directory
if [ -f output/iris.hdf5 ]; then
	mv output/output.hdf5 output.h5
fi
if [ -d output ]; then
	rm --recursive --force output
fi

# TODO: call the postprocessing script to generate the summary PDF
cd "<<root>>" || exit 1
module load anaconda3/2023.07-2
echo "python3 python/postprocess_iris_run <<name>> --status=${result}"

exit ${exit_code}
