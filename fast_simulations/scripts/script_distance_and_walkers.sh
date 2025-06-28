#!/bin/bash

# === OAR directives ===
#OAR -t night
#OAR -l walltime=12:00:00
#OAR -O OAR_%jobid%.out
#OAR -E OAR_%jobid%.err
#OAR -n distance-target

# === Parametri della simulazione ===
# Nome file CSV unico con PID e timestamp
SIM_CSV_FILE="distance_pid_$$-$(date +%s).csv"
SIM_TRIALS=20  # Numero di iterazioni per la simulazione

# === Nomi dei file ===
SOURCE_FILE="variable_distance_multi_target.c"
LIB_FILE="func.c"
EXECUTABLE_NAME="distance_$$"  # Eseguibile con nome unico

# === Caricamento modulo GCC ===
module load gcc  # Assicurati che questo modulo esista nel tuo cluster

# === Compilazione ===
echo "Compilazione di $SOURCE_FILE..."
gcc "$SOURCE_FILE" "$LIB_FILE" -o "$EXECUTABLE_NAME" -lm
if [ $? -ne 0 ]; then
    echo "Errore: Compilazione fallita per $SOURCE_FILE."
    exit 1
fi
echo "Compilazione riuscita. Eseguibile: $EXECUTABLE_NAME"

# === Esecuzione della simulazione ===
echo "Esecuzione della simulazione..."
"./$EXECUTABLE_NAME" "$SIM_CSV_FILE" "$SIM_TRIALS"
if [ $? -ne 0 ]; then
    echo "Errore: L'esecuzione della simulazione è fallita."
else
    echo "Simulazione completata con successo."
fi

# === Pulizia ===
echo "Cancellazione dell'eseguibile $EXECUTABLE_NAME..."
rm -f "./$EXECUTABLE_NAME"
if [ $? -ne 0 ]; then
    echo "Avviso: Impossibile cancellare l'eseguibile $EXECUTABLE_NAME."
fi

# === Fine ===
echo "Data e ora di fine: $(date)"
echo "--- Job OAR completato ---"
