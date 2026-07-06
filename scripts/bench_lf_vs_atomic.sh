#!/usr/bin/env bash
#
# Lock-free vs Atomic comparison with full ncu metrics.
#
# Reruns the lock-free (cell-parallel) vs atomic (edge-parallel) benchmark across
# the documented resolutions and, for every (resolution, strategy) pair, collects:
#   * wall-clock ns/cell/iter (benchmark's own averaged timer)
#   * per-kernel time via nsys   (x-sweep, y-sweep, wet/dry, initNewCells)
#   * per-kernel achieved occupancy, registers/thread and DRAM bandwidth via ncu
# into a single CSV that mirrors the block-size sweep produced by
# scripts/bench_blocksizes.sh.
#
# Configuration can be overridden via environment variables:
#
#   RESOLUTIONS   cell sizes to test              (default "2000m 1000m 500m")
#   STRATEGIES    sync strategies to test         (default "lock-free atomic")
#   RUNS          averaged runs for wall-clock    (default 10)
#   PROFILE_RUNS  runs performed under nsys       (default 2)
#   END_TIME      simulation end time in seconds  (default 10800)
#   BLOCK_SIZE    CUDA block width                (default 16)
#   BIN           benchmark binary                (default ./build/benchmark_cuda)
#   OUTDIR        output directory                (default data/lf_vs_atomic_sweep)
#   BUILD         set to 1 to build first         (default 0)
#   NCU           set to 0 to skip ncu profiling  (default 1)
#   NCU_BIN       ncu binary                      (default ncu)
#   NCU_END_TIME  sim end time for ncu pass       (default 60)
#   NCU_SKIP      kernel launches to skip         (default 6)
#   NCU_LAUNCHES  kernel launches to profile      (default 18)
#
set -euo pipefail

RESOLUTIONS="${RESOLUTIONS:-2000m 1000m 500m}"
STRATEGIES="${STRATEGIES:-lock-free atomic}"
RUNS="${RUNS:-10}"
PROFILE_RUNS="${PROFILE_RUNS:-2}"
END_TIME="${END_TIME:-10800}"
BLOCK_SIZE="${BLOCK_SIZE:-16}"
BIN="${BIN:-./build/benchmark_cuda}"
OUTDIR="${OUTDIR:-data/lf_vs_atomic_sweep}"
BUILD="${BUILD:-0}"
NCU="${NCU:-1}"
NCU_BIN="${NCU_BIN:-ncu}"
NCU_END_TIME="${NCU_END_TIME:-60}"
NCU_SKIP="${NCU_SKIP:-6}"
NCU_LAUNCHES="${NCU_LAUNCHES:-18}"

# occupancy + registers/thread + DRAM bandwidth metrics collected by ncu.
NCU_OCC_METRIC="sm__warps_active.avg.pct_of_peak_sustained_active"
NCU_REG_METRIC="launch__registers_per_thread"
NCU_BW_PCT_METRIC="dram__throughput.avg.pct_of_peak_sustained_elapsed"
NCU_BW_BPS_METRIC="dram__bytes.sum.per_second"

# --- sanity checks ----------------------------------------------------------
command -v nsys >/dev/null 2>&1 || {
  echo "error: nsys not found on PATH (are you inside the nix-shell?)" >&2
  exit 1
}

if [ "$NCU" = "1" ] && ! command -v "$NCU_BIN" >/dev/null 2>&1; then
  echo "warning: '$NCU_BIN' not found on PATH; disabling ncu profiling." >&2
  NCU=0
fi

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
CSV="$OUTDIR/lf_vs_atomic_results.csv"
echo "resolution,strategy,ns_per_cell_iter,total_kernel_ms_per_run,xsweep_ms_per_run,ysweep_ms_per_run,wetdry_ms_per_run,initcells_ms_per_run,xsweep_occ_pct,xsweep_regs,xsweep_dram_pct,xsweep_dram_gbps,ysweep_occ_pct,ysweep_regs,ysweep_dram_pct,ysweep_dram_gbps,wetdry_occ_pct,wetdry_regs,wetdry_dram_pct,wetdry_dram_gbps,initcells_occ_pct,initcells_regs,initcells_dram_pct,initcells_dram_gbps" > "$CSV"

# Extract "Total Time (ns)" (column 2 of the CSV report) for rows whose kernel
# name matches the given regex; sum them and divide by PROFILE_RUNS.
sum_kernel_ns() {
  local csv="$1" pattern="$2"
  echo "$csv" | awk -F',' -v pat="$pattern" '
    NR > 1 && $0 ~ pat { gsub(/"/, "", $2); s += $2 }
    END { printf "%.0f", s + 0 }'
}

# Average an ncu metric over all profiled launches of kernels whose name matches
# the given regex. Handles quoted CSV fields (gawk FPAT). Prints NA if missing.
ncu_metric() {
  local csv="$1" kpat="$2" mname="$3"
  echo "$csv" | awk -v FPAT='([^,]*)|("[^"]*")' -v kpat="$kpat" -v mname="$mname" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        g = $i; gsub(/"/, "", g)
        if (g == "Kernel Name")  kc = i
        if (g == "Metric Name")  mc = i
        if (g == "Metric Value") vc = i
      }
      next
    }
    kc && mc && vc {
      kn = $kc; gsub(/"/, "", kn)
      mn = $mc; gsub(/"/, "", mn)
      mv = $vc; gsub(/[",]/, "", mv)
      if (kn ~ kpat && mn == mname) { s += mv; c++ }
    }
    END { if (c > 0) printf "%.2f", s / c; else printf "NA" }'
}

for RES in $RESOLUTIONS; do
  for STRAT in $STRATEGIES; do
    TAG="${RES}_${STRAT}"
    echo
    echo "==================================================================="
    echo " resolution ${RES}   strategy ${STRAT}   (block ${BLOCK_SIZE}x${BLOCK_SIZE})"
    echo "==================================================================="

    # 1) wall-clock timing: the benchmark's own averaged timer.
    BENCH_LOG="$OUTDIR/bench_${TAG}.log"
    "$BIN" "$RUNS" "$END_TIME" \
        --resolution "$RES" --block-size "$BLOCK_SIZE" --sync-strategy "$STRAT" \
        | tee "$BENCH_LOG"

    NS_PER_CELL=$(awk '/^time per cell and iteration in ns:/ { v = $NF } END { print v }' "$BENCH_LOG")

    # 2) nsys profile (fewer runs; profiling adds overhead).
    REP="$OUTDIR/${TAG}"
    echo ">> profiling with nsys ($PROFILE_RUNS run(s)) ..."
    nsys profile --force-overwrite true -o "$REP" \
        "$BIN" "$PROFILE_RUNS" "$END_TIME" \
        --resolution "$RES" --block-size "$BLOCK_SIZE" --sync-strategy "$STRAT" \
        >/dev/null

    # 3) pull the per-kernel totals out of the nsys report as CSV.
    KERN_CSV=$(nsys stats --force-export=true --report cuda_gpu_kern_sum --format csv "${REP}.nsys-rep" 2>/dev/null || true)

    if [ -z "$KERN_CSV" ]; then
      echo "warning: could not read cuda_gpu_kern_sum for $TAG" >&2
      echo "${RES},${STRAT},${NS_PER_CELL:-NA},NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> "$CSV"
      continue
    fi

    TOTAL_NS=$(sum_kernel_ns "$KERN_CSV" '::cuda::')
    XS_NS=$(sum_kernel_ns "$KERN_CSV" 'computeXSweep')
    YS_NS=$(sum_kernel_ns "$KERN_CSV" 'computeYSweep')
    WD_NS=$(sum_kernel_ns "$KERN_CSV" 'applyWetDryThreshold')
    IC_NS=$(sum_kernel_ns "$KERN_CSV" 'initNewCells')

    # ns over all profiled runs -> ms per single run.
    to_ms_per_run() { awk -v v="$1" -v r="$PROFILE_RUNS" 'BEGIN { printf "%.3f", (v / r) / 1e6 }'; }
    TOTAL_MS=$(to_ms_per_run "$TOTAL_NS")
    XS_MS=$(to_ms_per_run "$XS_NS")
    YS_MS=$(to_ms_per_run "$YS_NS")
    WD_MS=$(to_ms_per_run "$WD_NS")
    IC_MS=$(to_ms_per_run "$IC_NS")

    # 4) ncu profile: achieved occupancy + registers/thread + DRAM bandwidth per kernel.
    XS_OCC=NA; XS_REG=NA; XS_BWP=NA; XS_BWG=NA
    YS_OCC=NA; YS_REG=NA; YS_BWP=NA; YS_BWG=NA
    WD_OCC=NA; WD_REG=NA; WD_BWP=NA; WD_BWG=NA
    IC_OCC=NA; IC_REG=NA; IC_BWP=NA; IC_BWG=NA
    if [ "$NCU" = "1" ]; then
      echo ">> profiling with ncu (skip $NCU_SKIP, count $NCU_LAUNCHES launches) ..."
      NCU_OUT=$("$NCU_BIN" --csv --target-processes all \
          --metrics "${NCU_OCC_METRIC},${NCU_REG_METRIC},${NCU_BW_PCT_METRIC},${NCU_BW_BPS_METRIC}" \
          --launch-skip "$NCU_SKIP" --launch-count "$NCU_LAUNCHES" \
          "$BIN" 1 "$NCU_END_TIME" \
          --resolution "$RES" --block-size "$BLOCK_SIZE" --sync-strategy "$STRAT" \
          2>/dev/null || true)

      if [ -z "$NCU_OUT" ]; then
        echo "warning: ncu produced no output for $TAG" >&2
      else
        # strip everything before the CSV header ncu prints ("ID",...).
        NCU_CSV=$(echo "$NCU_OUT" | sed -n '/"ID"/,$p')

        # convert an averaged bytes/s figure to GB/s (1e9 bytes); passes NA through.
        to_gbps() { awk -v v="$1" 'BEGIN { if (v == "NA") print "NA"; else printf "%.2f", v / 1e9 }'; }

        XS_OCC=$(ncu_metric "$NCU_CSV" 'computeXSweep'        "$NCU_OCC_METRIC")
        XS_REG=$(ncu_metric "$NCU_CSV" 'computeXSweep'        "$NCU_REG_METRIC")
        XS_BWP=$(ncu_metric "$NCU_CSV" 'computeXSweep'        "$NCU_BW_PCT_METRIC")
        XS_BWG=$(to_gbps "$(ncu_metric "$NCU_CSV" 'computeXSweep'        "$NCU_BW_BPS_METRIC")")
        YS_OCC=$(ncu_metric "$NCU_CSV" 'computeYSweep'        "$NCU_OCC_METRIC")
        YS_REG=$(ncu_metric "$NCU_CSV" 'computeYSweep'        "$NCU_REG_METRIC")
        YS_BWP=$(ncu_metric "$NCU_CSV" 'computeYSweep'        "$NCU_BW_PCT_METRIC")
        YS_BWG=$(to_gbps "$(ncu_metric "$NCU_CSV" 'computeYSweep'        "$NCU_BW_BPS_METRIC")")
        WD_OCC=$(ncu_metric "$NCU_CSV" 'applyWetDryThreshold' "$NCU_OCC_METRIC")
        WD_REG=$(ncu_metric "$NCU_CSV" 'applyWetDryThreshold' "$NCU_REG_METRIC")
        WD_BWP=$(ncu_metric "$NCU_CSV" 'applyWetDryThreshold' "$NCU_BW_PCT_METRIC")
        WD_BWG=$(to_gbps "$(ncu_metric "$NCU_CSV" 'applyWetDryThreshold' "$NCU_BW_BPS_METRIC")")
        IC_OCC=$(ncu_metric "$NCU_CSV" 'initNewCells'         "$NCU_OCC_METRIC")
        IC_REG=$(ncu_metric "$NCU_CSV" 'initNewCells'         "$NCU_REG_METRIC")
        IC_BWP=$(ncu_metric "$NCU_CSV" 'initNewCells'         "$NCU_BW_PCT_METRIC")
        IC_BWG=$(to_gbps "$(ncu_metric "$NCU_CSV" 'initNewCells'         "$NCU_BW_BPS_METRIC")")
      fi
    fi

    echo "${RES},${STRAT},${NS_PER_CELL:-NA},${TOTAL_MS},${XS_MS},${YS_MS},${WD_MS},${IC_MS},${XS_OCC},${XS_REG},${XS_BWP},${XS_BWG},${YS_OCC},${YS_REG},${YS_BWP},${YS_BWG},${WD_OCC},${WD_REG},${WD_BWP},${WD_BWG},${IC_OCC},${IC_REG},${IC_BWP},${IC_BWG}" >> "$CSV"
  done
done

echo
echo "=== results: $CSV ==="
if command -v column >/dev/null 2>&1; then
  column -s, -t "$CSV"
else
  cat "$CSV"
fi

echo
if [ "$NCU" = "1" ]; then
  echo "note: per-kernel ncu columns, averaged over the profiled launches:"
  echo "      *_occ_pct   = achieved occupancy   (${NCU_OCC_METRIC})"
  echo "      *_regs      = registers/thread     (${NCU_REG_METRIC})"
  echo "      *_dram_pct  = DRAM BW % of peak    (${NCU_BW_PCT_METRIC})"
  echo "      *_dram_gbps = DRAM bandwidth GB/s  (${NCU_BW_BPS_METRIC})"
  echo "      initcells_* is atomic-only (the lock-free solver has no initNewCells"
  echo "      kernel), so those columns are NA for lock-free rows."
  echo "      DRAM metrics may require multiple kernel replays; lower NCU_LAUNCHES"
  echo "      if the ncu pass is too slow. Set NCU=0 to skip the ncu pass."
else
  echo "note: ncu profiling was skipped (NCU=0 or ncu unavailable); the occupancy,"
  echo "      registers/thread and DRAM bandwidth columns are NA. Re-run with NCU=1"
  echo "      and ncu on PATH."
fi
