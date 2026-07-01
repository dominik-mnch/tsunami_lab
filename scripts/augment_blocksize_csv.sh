#!/usr/bin/env bash
#
# Task 4 - augment the block-size sweep CSV from existing nsys reports.
#
# Reads an existing blocksize_results.csv (produced by bench_blocksizes.sh) and,
# for each block size, pulls extra metrics out of the matching bs<N>.nsys-rep
# file that is already on disk - no re-profiling required. The extra columns are
# appended and written to a new *_augmented.csv (the original is left untouched).
#
# Extracted per report:
#   h2d_ms      host-to-device copy time (ms, summed over the profiled runs)
#   h2d_mb      host-to-device volume (MB, summed)
#   h2d_gbps    achieved H2D transfer bandwidth (GB/s) = h2d_mb*1e6 / h2d_ns
#   sync_pct    share of CUDA API time in cudaDeviceSynchronize
#   launch_pct  share of CUDA API time in cudaLaunchKernel
#
# Overridable via environment variables:
#   OUTDIR  sweep directory        (default data/blocksize_sweep)
#   CSV     input CSV              (default $OUTDIR/blocksize_results.csv)
#   OUT     output CSV             (default ${CSV%.csv}_augmented.csv)
#
set -euo pipefail

OUTDIR="${OUTDIR:-data/blocksize_sweep}"
CSV="${CSV:-$OUTDIR/blocksize_results.csv}"
OUT="${OUT:-${CSV%.csv}_augmented.csv}"

command -v nsys >/dev/null 2>&1 || {
  echo "error: nsys not found on PATH (are you inside the nix-shell?)" >&2
  exit 1
}
[ -f "$CSV" ] || {
  echo "error: input CSV '$CSV' not found (run bench_blocksizes.sh first)" >&2
  exit 1
}

# --- helpers: pull single aggregate values out of an nsys report ------------

# Sum of host-to-device copy time (ns).
h2d_ns() {
  nsys stats --report cuda_gpu_mem_time_sum --format csv "$1" 2>/dev/null \
    | awk -F',' '/Host-to-Device/ { gsub(/"/, "", $2); s += $2 } END { printf "%.0f", s + 0 }'
}

# Sum of host-to-device volume (MB).
h2d_mb() {
  nsys stats --report cuda_gpu_mem_size_sum --format csv "$1" 2>/dev/null \
    | awk -F',' '/Host-to-Device/ { gsub(/"/, "", $1); s += $1 } END { printf "%.3f", s + 0 }'
}

# Time (%) column for the API call whose name matches the given pattern.
api_pct() {
  nsys stats --report cuda_api_sum --format csv "$1" 2>/dev/null \
    | awk -F',' -v pat="$2" '$0 ~ pat { gsub(/"/, "", $1); v = $1 } END { printf "%.1f", v + 0 }'
}

# --- join the extra columns onto every row of the existing CSV --------------

NEW_HEADER="h2d_ms,h2d_mb,h2d_gbps,sync_pct,launch_pct"

{
  read -r header
  echo "${header},${NEW_HEADER}"

  while IFS= read -r line || [ -n "$line" ]; do
    [ -z "$line" ] && continue

    N="${line%%,*}"                       # block_size is the first CSV field
    REP="$OUTDIR/bs${N}.nsys-rep"

    if [ ! -f "$REP" ]; then
      echo "warning: report '$REP' missing; leaving row $N blank" >&2
      echo "${line},NA,NA,NA,NA,NA"
      continue
    fi

    NS=$(h2d_ns "$REP")
    MB=$(h2d_mb "$REP")
    MS=$(awk -v v="$NS" 'BEGIN { printf "%.3f", v / 1e6 }')
    GBPS=$(awk -v mb="$MB" -v ns="$NS" 'BEGIN { if (ns > 0) printf "%.1f", mb * 1e6 / ns; else printf "NA" }')
    SYNC=$(api_pct "$REP" 'cudaDeviceSynchronize')
    LAUNCH=$(api_pct "$REP" 'cudaLaunchKernel')

    echo "${line},${MS},${MB},${GBPS},${SYNC},${LAUNCH}"
  done
} < "$CSV" > "$OUT"

echo "=== augmented: $OUT ==="
if command -v column >/dev/null 2>&1; then
  column -s, -t "$OUT"
else
  cat "$OUT"
fi

echo
echo "note: h2d_gbps is the host-to-device *transfer* bandwidth and is expected to"
echo "      be roughly constant across block sizes (the copies do not depend on"
echo "      block size). For kernel DRAM throughput you still need nsys"
echo "      --gpu-metrics-device or ncu. The sync_pct / launch_pct columns do vary"
echo "      with block size and show how launch overhead shrinks as blocks grow."
