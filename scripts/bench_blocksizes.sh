#!/usr/bin/env bash
#
# Task 4 - Thread Count Benchmarking.
#
# Sweeps the CUDA block size, runs the GPU benchmark for each, profiles it with
# nsys, and extracts the relevant kernel timings into a single CSV that can be
# pasted straight into the Task 4 comparison table.
#
# Configuration can be overridden via environment variables:
#
#   BLOCK_SIZES   block widths to test           (default "1 2 4 8 16 32")
#   RUNS          averaged runs for wall-clock    (default 10)
#   PROFILE_RUNS  runs performed under nsys       (default 2)
#   END_TIME      simulation end time in seconds  (default 10800)
#   RESOLUTION    cell size in metres             (default 1000m)
#   SYNC          sync strategy lock-free|atomic  (default lock-free)
#   BIN           benchmark binary                (default ./build/benchmark_cuda)
#   OUTDIR        output directory                (default data/blocksize_sweep)
#   BUILD         set to 1 to build first         (default 0)
#
set -euo pipefail

BLOCK_SIZES="${BLOCK_SIZES:-1 2 4 8 16 32}"
RUNS="${RUNS:-10}"
PROFILE_RUNS="${PROFILE_RUNS:-2}"
END_TIME="${END_TIME:-10800}"
RESOLUTION="${RESOLUTION:-1000m}"
SYNC="${SYNC:-lock-free}"
BIN="${BIN:-./build/benchmark_cuda}"
OUTDIR="${OUTDIR:-data/blocksize_sweep}"
BUILD="${BUILD:-0}"

# --- sanity checks ----------------------------------------------------------
command -v nsys >/dev/null 2>&1 || {
  echo "error: nsys not found on PATH (are you inside the nix-shell?)" >&2
  exit 1
}

if [ "$BUILD" = "1" ]; then
  echo ">> building (scons mode=release) ..."
  scons mode=release -j4
fi

[ -x "$BIN" ] || {
  echo "error: benchmark binary '$BIN' not found or not executable." >&2
  echo "       build it first, e.g. BUILD=1 $0" >&2
  exit 1
}

mkdir -p "$OUTDIR"
CSV="$OUTDIR/blocksize_results.csv"
echo "block_size,threads_per_block,ns_per_cell_iter,total_kernel_ms_per_run,xsweep_ms_per_run,ysweep_ms_per_run,wetdry_ms_per_run" > "$CSV"

# Extract "Total Time (ns)" (column 2 of the CSV report) for rows whose kernel
# name matches the given regex; sum them and divide by PROFILE_RUNS.
sum_kernel_ns() {
  local csv="$1" pattern="$2"
  echo "$csv" | awk -F',' -v pat="$pattern" '
    NR > 1 && $0 ~ pat { gsub(/"/, "", $2); s += $2 }
    END { printf "%.0f", s + 0 }'
}

for N in $BLOCK_SIZES; do
  TPB=$(( N * N ))
  echo
  echo "==================================================================="
  echo " block size ${N}  (${N}x${N} = ${TPB} threads/block)"
  echo "==================================================================="

  # 1) wall-clock timing: the benchmark's own averaged timer.
  BENCH_LOG="$OUTDIR/bench_bs${N}.log"
  "$BIN" "$RUNS" "$END_TIME" \
      --resolution "$RESOLUTION" --block-size "$N" --sync-strategy "$SYNC" \
      | tee "$BENCH_LOG"

  NS_PER_CELL=$(awk '/^time per cell and iteration in ns:/ { v = $NF } END { print v }' "$BENCH_LOG")

  # 2) nsys profile (fewer runs; profiling adds overhead).
  REP="$OUTDIR/bs${N}"
  echo ">> profiling with nsys ($PROFILE_RUNS run(s)) ..."
  nsys profile --force-overwrite true -o "$REP" \
      "$BIN" "$PROFILE_RUNS" "$END_TIME" \
      --resolution "$RESOLUTION" --block-size "$N" --sync-strategy "$SYNC" \
      >/dev/null

  # 3) pull the per-kernel totals out of the nsys report as CSV.
  KERN_CSV=$(nsys stats --report cuda_gpu_kern_sum --format csv "${REP}.nsys-rep" 2>/dev/null || true)

  if [ -z "$KERN_CSV" ]; then
    echo "warning: could not read cuda_gpu_kern_sum for block size $N" >&2
    echo "${N},${TPB},${NS_PER_CELL:-NA},NA,NA,NA,NA" >> "$CSV"
    continue
  fi

  TOTAL_NS=$(sum_kernel_ns "$KERN_CSV" '::cuda::')
  XS_NS=$(sum_kernel_ns "$KERN_CSV" 'computeXSweep')
  YS_NS=$(sum_kernel_ns "$KERN_CSV" 'computeYSweep')
  WD_NS=$(sum_kernel_ns "$KERN_CSV" 'applyWetDryThreshold')

  # ns over all profiled runs -> ms per single run.
  to_ms_per_run() { awk -v v="$1" -v r="$PROFILE_RUNS" 'BEGIN { printf "%.3f", (v / r) / 1e6 }'; }
  TOTAL_MS=$(to_ms_per_run "$TOTAL_NS")
  XS_MS=$(to_ms_per_run "$XS_NS")
  YS_MS=$(to_ms_per_run "$YS_NS")
  WD_MS=$(to_ms_per_run "$WD_NS")

  echo "${N},${TPB},${NS_PER_CELL:-NA},${TOTAL_MS},${XS_MS},${YS_MS},${WD_MS}" >> "$CSV"
done

echo
echo "=== results: $CSV ==="
if command -v column >/dev/null 2>&1; then
  column -s, -t "$CSV"
else
  cat "$CSV"
fi

echo
echo "note: occupancy and registers/thread are NOT available from nsys stats."
echo "      registers/thread: rebuild with nvcc '-Xptxas -v' and read 'Used N registers'."
echo "      occupancy:        use 'ncu --metrics sm__warps_active.avg.pct_of_peak_sustained_active'"
echo "                        (falls back to nsys --gpu-metrics if ncu is blocked on WDDM)."
