#!/bin/bash

#OAR -t night
#OAR -l walltime=12:00:00
#OAR -p cluster=musa
#OAR -O OAR_%jobid%.out
#OAR -E OAR_%jobid%.err
#OAR -n 3d-levy-walks

# --- Display resource info ---
hostname
echo "Starting job on $(hostname) at $(date)"
echo "Current working directory: $(pwd)"

# --- Activate Conda Environment ---
# 1. Carica il modulo conda (questa è la tua prima riga)
module load conda

# 2. Verifica che 'conda' sia ora nel PATH (debugging)
echo "Conda path after module load: $(which conda)"
echo "Conda version after module load: $(conda --version)"

# 3. Attiva l'ambiente desiderato
conda activate levykernel

# 4. Check if activation was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment 'levykernel'. Exiting."
    exit 1
fi

echo "Conda environment 'levykernel' activated."
echo "Python interpreter used: $(which python)"
echo "Python version: $(python --version)"

# --- Execute the Python script ---
PYTHON_SCRIPT="3d-experiments.py"
FULL_PYTHON_SCRIPT_PATH="$(pwd)/${PYTHON_SCRIPT}"

if [ -f "${FULL_PYTHON_SCRIPT_PATH}" ]; then
    echo "Executing Python script: ${FULL_PYTHON_SCRIPT_PATH}"
    python "${FULL_PYTHON_SCRIPT_PATH}"
    if [ $? -ne 0 ]; then
        echo "Error: Python script '${PYTHON_SCRIPT}' failed. Exiting with error."
        exit 1
    else
        echo "Python script '${PYTHON_SCRIPT}' completed successfully."
    fi
else
    echo "Error: Python script '${PYTHON_SCRIPT}' not found at '${FULL_PYTHON_SCRIPT_PATH}'. Exiting."
    exit 1
fi

echo "Job finished at $(date)"