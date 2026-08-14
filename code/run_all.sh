#!/usr/bin/env bash
# run_all.sh
# Master script: run ALL simulations then ALL table/figure scripts
#
#   chmod +x run_all.sh
#   nohup ./run_all.sh > run_all_master.log 2>&1 &

set -uo pipefail

MAX_JOBS=5           # max concurrent simulation jobs
MAX_CORES=102        # max CPU cores per Rscript job (each job uses up to MAX_CORES workers)
POLL_INTERVAL=30     # seconds between status checks

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LOG_DIR="$SCRIPT_DIR/../output/logs"
mkdir -p "$LOG_DIR"
export N_CORES=$MAX_CORES   # propagate core limit to all Rscript subprocesses

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
MASTER_LOG="$LOG_DIR/master_${TIMESTAMP}.log"

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg" | tee -a "$MASTER_LOG"
}

log "=========================================="
log "Starting job scheduler -- MAX_JOBS=$MAX_JOBS, MAX_CORES=$MAX_CORES"
log "Working directory: $SCRIPT_DIR"
log "Log directory: $LOG_DIR"
log "=========================================="

# Check for PPMI data files
for f in Mimic_PPMI_data.RData Fitting_PPMI.RData; do
    if [ -f "$f" ]; then
        log "Data file found: $f"
    else
        log "WARNING: $f not found (MimicPPMI tasks may fail)"
    fi
done

# ========================== Build job list ==========================
JOBS=()

# Simulation_Main.R: params 1-8
for i in $(seq 1 8); do
    JOBS+=("Rscript Simulation_Main.R $i")
done

# Simulation_Addition.R: params 1-16
for i in $(seq 1 16); do
    JOBS+=("Rscript Simulation_Addition.R $i")
done

# Simulation_Suggestion.R: params 1-8
for i in $(seq 1 8); do
    JOBS+=("Rscript Simulation_Suggestion.R $i")
done

# Simulation_Diagnosis.R: no params
JOBS+=("Rscript Simulation_Diagnosis.R")

# Simulation_ModelSelection.R: params 1-6
for i in $(seq 1 6); do
    JOBS+=("Rscript Simulation_ModelSelection.R $i")
done

# Simulation_MimicPPMI.R: params 3-7
for i in $(seq 3 7); do
    JOBS+=("Rscript Simulation_MimicPPMI.R $i")
done

# Simulation_MimicPPMI_ModelSelection.R: no params
JOBS+=("Rscript Simulation_MimicPPMI_ModelSelection.R")

TOTAL_JOBS=${#JOBS[@]}
log "Total simulation jobs: $TOTAL_JOBS"

# ========================== Job pool management ==========================
declare -A PID_TO_CMD
declare -A PID_TO_LOG
declare -A PID_TO_START
JOB_INDEX=0
COMPLETED=0
FAILED=0
FAILED_LIST=()

# Submit next job
submit_next() {
    if [ $JOB_INDEX -lt $TOTAL_JOBS ]; then
        local cmd="${JOBS[$JOB_INDEX]}"
        local job_num=$((JOB_INDEX + 1))
        local safe_name
        safe_name=$(echo "$cmd" | sed 's/Rscript //; s/ /_/g; s/%/pct/g')
        local job_log="$LOG_DIR/${safe_name}_${TIMESTAMP}.log"

        log "[$job_num/$TOTAL_JOBS] Starting: $cmd"

        eval "$cmd" > "$job_log" 2>&1 &
        local pid=$!

        PID_TO_CMD[$pid]="$cmd"
        PID_TO_LOG[$pid]="$job_log"
        PID_TO_START[$pid]=$(date +%s)

        JOB_INDEX=$((JOB_INDEX + 1))
    fi
}

format_duration() {
    local seconds=$1
    local hours=$((seconds / 3600))
    local minutes=$(( (seconds % 3600) / 60 ))
    local secs=$((seconds % 60))
    printf "%dh%02dm%02ds" $hours $minutes $secs
}

# ========================== Main loop ==========================

# Fill initial job pool
for _ in $(seq 1 $MAX_JOBS); do
    submit_next
done

# Monitor and refill
while [ ${#PID_TO_CMD[@]} -gt 0 ]; do
    sleep "$POLL_INTERVAL"

    for pid in "${!PID_TO_CMD[@]}"; do
        if ! kill -0 "$pid" 2>/dev/null; then
            wait "$pid" 2>/dev/null
            exit_code=$?

            local_cmd="${PID_TO_CMD[$pid]}"
            local_log="${PID_TO_LOG[$pid]}"
            local_start="${PID_TO_START[$pid]}"
            local_end=$(date +%s)
            local_dur=$(format_duration $((local_end - local_start)))

            unset "PID_TO_CMD[$pid]"
            unset "PID_TO_LOG[$pid]"
            unset "PID_TO_START[$pid]"

            COMPLETED=$((COMPLETED + 1))

            if [ "$exit_code" -eq 0 ]; then
                log "[Done $COMPLETED/$TOTAL_JOBS] OK: $local_cmd (time: $local_dur)"
            else
                FAILED=$((FAILED + 1))
                FAILED_LIST+=("$local_cmd")
                log "[Fail $COMPLETED/$TOTAL_JOBS] FAILED: $local_cmd (exit $exit_code, time: $local_dur)"
                log "  Log: $local_log"
            fi

            submit_next
        fi
    done

    running=${#PID_TO_CMD[@]}
    remaining=$((TOTAL_JOBS - JOB_INDEX))
    log "Status: running=$running, completed=$COMPLETED/$TOTAL_JOBS, queued=$remaining, failed=$FAILED"
done

log "=========================================="
log "All simulation jobs done. Success=$((COMPLETED - FAILED)), Failed=$FAILED"
if [ $FAILED -gt 0 ]; then
    log "Failed jobs:"
    for cmd in "${FAILED_LIST[@]}"; do
        log "  - $cmd"
    done
fi
log "=========================================="

# ========================== Phase 2: Table scripts ==========================
log ""
log "=========================================="
log "Running Table scripts (sequential)..."
log "=========================================="

TABLE_SCRIPTS=(
    TableS1_MainSimulation.R
    TableS2_ModelSelection.R
    TableS3_MimicPPMI.R
    TableS4_AdditionalRobustness.R
    TableS5_EstimationStrategies.R
    MimicPPMI_ModelSelection_Summary.R
)
TABLE_FAILED=0

for script in "${TABLE_SCRIPTS[@]}"; do
    if [ ! -f "$script" ]; then
        log "Skipping $script (file not found)"
        continue
    fi

    local_log="$LOG_DIR/${script%.R}_${TIMESTAMP}.log"
    log "Running: Rscript $script"
    start_time=$(date +%s)

    if Rscript "$script" > "$local_log" 2>&1; then
        end_time=$(date +%s)
        dur=$(format_duration $((end_time - start_time)))
        log "  OK: $script done (time: $dur)"
    else
        end_time=$(date +%s)
        dur=$(format_duration $((end_time - start_time)))
        TABLE_FAILED=$((TABLE_FAILED + 1))
        log "  FAILED: $script (time: $dur), log: $local_log"
    fi
done

log "=========================================="
log "All done!"
log "  Simulations: success=$((COMPLETED - FAILED))/$TOTAL_JOBS, failed=$FAILED"
log "  Tables: done=$((${#TABLE_SCRIPTS[@]} - TABLE_FAILED))/${#TABLE_SCRIPTS[@]}, failed=$TABLE_FAILED"
log "  All logs saved in: $LOG_DIR/"
log "=========================================="
