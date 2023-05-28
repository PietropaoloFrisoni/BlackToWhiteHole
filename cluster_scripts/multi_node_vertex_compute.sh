#!/bin/bash
#SBATCH -A def-vidotto
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=3-0:00:00
#SBATCH --job-name=multi_node_vertex_compute
#SBATCH --output=multi_node_vertex_compute.log
#SBATCH --error=multi_node_vertex_compute.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=pfrisoni@uwo.ca

echo "Running on: $SLURM_NODELIST"
echo

# start commands

BASEDIR=/home/frisus95/projects/def-vidotto/frisus95/sl2cfoam-next
DATADIR=${BASEDIR}/data_sl2cfoam

export LD_LIBRARY_PATH="${BASEDIR}/lib":$LD_LIBRARY_PATH

# number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

IMMIRZI=1.0
SHELLmin=0
SHELLmax=8

# Define the spin values to loop over
# (chosen taking into account triangular inequalities)
TSPINVALUES=(11 5 12 6 13 5 14 6 15 7 16 6 17 7 18 8 19 7 20 8)

# Loop over the values using two variables
for ((i = 0; i < ${#TSPINVALUES[@]}; i += 2)); do

      TJ0=${TSPINVALUES[i]}
      TJPM=${TSPINVALUES[i + 1]}

      for ((CURRENTSHELL = $SHELLmin; CURRENTSHELL <= $SHELLmax; CURRENTSHELL += 1)); do

            now=$(date)
            echo
            echo "Starting Lorentzian fulltensor [ IMMIRZI = ${IMMIRZI} TJ0 = ${TJ0} TJPM = ${TJPM}, shells = ${CURRENTSHELL} ]... (now: $now)"

            echo "TJ0: $TJ0, TJPM: $TJPM"

            $BASEDIR/bin/vertex-fulltensor -V -h -m 2000 $DATADIR $IMMIRZI $TJ0,$TJ0,$TJ0,$TJ0,$TJPM,$TJPM,$TJPM,$TJPM,$TJPM,$TJPM $CURRENTSHELL

            now=$(date)
            echo "... done (now: $now)"
            echo

      done # shell cycle

done # boundary spin cycle

echo
echo "All completed."
