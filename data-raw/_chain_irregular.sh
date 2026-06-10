#!/usr/bin/env bash
#
# Chain the three irregular-grid benchmarks: wait for adlift to finish
# (already running), then run nlt, then run rLifting irregular. Each step
# polls for its predecessor's .rda output before launching, so the chain
# survives a sleep/wake cycle but NOT a machine reboot.
#
# Launch detached:
#   cd /home/augusto/Documentos/projetos/rLifting
#   nohup bash data-raw/_chain_irregular.sh > /tmp/chain_irr.log 2>&1 &
#   disown
#
# Monitor:
#   tail -f /tmp/chain_irr.log
#
# Abort:
#   pkill -f _chain_irregular.sh
#   pkill -f generate_nlt_irregular
#   pkill -f generate_rlifting_irregular

set -u

PROJECT_DIR="/home/augusto/Documentos/projetos/rLifting"
cd "$PROJECT_DIR" || exit 1

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

wait_for_file() {
  local f="$1"
  log "Waiting for $f to appear (polling every 5 min)..."
  while [ ! -f "$f" ]; do
    sleep 300
  done
  log "Detected $f (size: $(stat -c%s "$f") bytes)"
}

run_bench() {
  local script="$1"
  local label="$2"
  log "Starting $label benchmark: Rscript $script"
  Rscript "$script"
  local rc=$?
  if [ $rc -ne 0 ]; then
    log "ERROR: $label benchmark exited with code $rc; aborting chain."
    exit $rc
  fi
  log "$label benchmark finished cleanly."
}

log "Chain started (PID $$). Watching for adlift output."

# Step 1: wait for adlift (already running independently)
wait_for_file "data/benchmark_adlift_irregular.rda"

# Step 2: nlt
run_bench "data-raw/generate_nlt_irregular_benchmark.R" "nlt"
wait_for_file "data/benchmark_nlt_irregular.rda"

# Step 3: rlifting irregular
run_bench "data-raw/generate_rlifting_irregular_benchmark.R" "rlifting-irregular"
wait_for_file "data/benchmark_rlifting_irregular.rda"

log "All three irregular benchmarks complete. Chain exiting normally."
