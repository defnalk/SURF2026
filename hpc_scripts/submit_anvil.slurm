#!/bin/bash
#SBATCH --job-name=surf_lmp
#SBATCH --account=YOUR_ALLOCATION           # <-- ACCESS/Anvil allocation name
#SBATCH --partition=wholenode                # Anvil CPU partitions: shared | wholenode | wide
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128                # Anvil wholenode = 128 cores/node
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@caltech.edu

# -------------------------------------------------------------------
#  SURF 2026 — LAMMPS production run on Purdue Anvil
#  Usage:  sbatch submit_anvil.slurm <input_file.in>
#  Example: sbatch submit_anvil.slurm ../simulations/monovalent/prod_298K.in
# -------------------------------------------------------------------

set -euo pipefail

INPUT=${1:-in.lammps}
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: LAMMPS input file '$INPUT' not found." >&2
    exit 1
fi

# ---- Environment ---------------------------------------------------
module purge
module load gcc/11.2.0
module load openmpi/4.1.6
module load lammps/2Aug2023                 # adjust to installed version
module list                                 # record versions in slurm log

# ---- Run info ------------------------------------------------------
echo "Job ID        : $SLURM_JOB_ID"
echo "Job name      : $SLURM_JOB_NAME"
echo "Nodes         : $SLURM_JOB_NUM_NODES"
echo "Tasks         : $SLURM_NTASKS"
echo "Host          : $(hostname)"
echo "Submit dir    : $SLURM_SUBMIT_DIR"
echo "Input deck    : $INPUT"
echo "Started at    : $(date)"

cd "$SLURM_SUBMIT_DIR"

# ---- Launch LAMMPS -------------------------------------------------
# -sf omp/intel/kk can be enabled if the Anvil build supports it.
srun --mpi=pmi2 -n "$SLURM_NTASKS" lmp -in "$INPUT" -log log.lammps.$SLURM_JOB_ID

echo "Finished at   : $(date)"
