#!/bin/bash
#SBATCH --job-name=LILAC_<<name>>
#SBATCH --time=02:00:00
#SBATCH --partition=blizzard
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --qos=b-standard
#SBATCH --output=<<directory>>/lilac_%A.log
#SBATCH --mail-type=END

# Load LILAC and clean up the environment
module purge
module load lilac

# move to the directory where runs happen
cd "<<directory>>" || exit 1

echo "LILAC run '<<name>>' starts."
echo "$(date +'%m-%d %H:%M') | LILAC run '<<name>>' starts." >> ../../../runs.log

# run LILAC
lilac
exit_code=$?

# figure out if the run was successful, and if not, why not
if [ $exit_code -eq 0 ]; then
	if grep -q "negative temperature detected" "lilac_$SLURM_JOB_ID.log"; then
		exit_code=4
		result="failed"
		exit_string="LILAC run '<<name>>' fails with a negative temperature error."
	elif grep -q "run summary" "lilac_$SLURM_JOB_ID.log"; then
		exit_code=0
		result="completed"
		exit_string="LILAC run '<<name>>' exits successfully."
	else
		exit_code=3
		result="failed"
		exit_string="LILAC run '<<name>>' quits prematurely (see '<<directory>>/lilac_$SLURM_JOB_ID.log')."
	fi
else
	result="failed"
	exit_string="LILAC run '<<name>>' fails with error_code ${exit_code} (see '<<directory>>/lilac_$SLURM_JOB_ID.log')."
fi

# log the result to the slurm log and to runs.log
echo "${exit_string}"
echo "$(date +'%m-%d %H:%M') | ${exit_string}" >> ../../../runs.log

# Move files out of the "out" directory
if [ -f out/fort.13 ]; then
	mv "out/fort.13" "output.lpf"
fi
if [ -f out/lilac.hdf5 ]; then
	mv "out/lilac.hdf5" "output.h5"
fi
if [ -f out/lilac_output.txt ]; then
	mv "out/lilac_output.txt" "output.txt"
fi
if [ -d out ]; then
	rm --recursive --force out
fi

# call the postprocessing script to update the run_outputs.csv table and
cd "../../.."
module load anaconda3/2023.07-2
PYTHONPATH=. python3 python/postprocess_lilac_run.py "<<name>>" --status="${result}"

exit ${exit_code}
