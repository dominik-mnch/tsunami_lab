============================================================
CUDA Memory Layout Transformations: Column-Major and Tiled
============================================================

Memory Layout Frameworks
========================

To maximize global memory throughput on the GPU, data layout must align with hardware 
coalescing behaviors. Below are the three memory access models we evaluated:

Row-Major (Baseline)
--------------------
Consecutive threads in a warp step through adjacent column indices. Memory stride is 
governed by the total width of the grid, including ghost boundaries:

.. code-block:: c

   #define ROW_MAJOR(row, col, stride) ((row) * (stride) + (col))

Column-Major
------------
Consecutive elements in a column are contiguous in memory. This shifts the unit-stride 
to the vertical dimension ($y$-axis), optimizing sweeps that process data column-wise:

.. code-block:: c

   #define COLUMN_MAJOR(row, col, height) ((col) * (height) + (row))

Tiled 32x32 Assembly
--------------------
Data is broken down into a two-dimensional grid of independent $32 \times 32$ tiles. 
Within each tile, memory is row-major. This configuration localizes thread blocks to 
structured hardware boundaries:

.. code-block:: c

   const t_idx l_tileX = i_col / i_tileSize;
   const t_idx l_tileY = i_row / i_tileSize;
   const t_idx l_localX = i_col % i_tileSize;
   const t_idx l_localY = i_row % i_tileSize;
   const t_idx l_tileOffset = ((l_tileY * i_nTilesX + l_tileX) * i_tileSize * i_tileSize);
   return l_tileOffset + (l_localY * i_tileSize + l_localX);

---

Kernel Modifications & Dimensional Shifts
=========================================

1. Ghost Cell Boundary Enforcements
-----------------------------------

Boundary kernels are highly sensitive to layout orientation due to memory coalescing.

* **Left/Right Boundaries:**
  
  * *Column-Major:* Threads process a continuous vertical column vector. Consecutive thread IDs read and write to sequential addresses, creating perfectly coalesced global memory operations.
  * *Tiled:* Linear index mappings split these boundaries into individual segments across multiple coordinate blocks, requiring explicit evaluation via ``tileLinearIndex``.

* **Bottom/Top Boundaries:**
  
  * *Column-Major:* Because stride jumps across entire column heights (``l_col * i_height``), a linear row access patterns requires strided steps, impacting vectorization efficiency relative to Row-Major.

2. Finite Volume Sweeps (X and Y Passes)
----------------------------------------

The operator-splitting scheme introduces contrasting performance footprints between layouts:

.. list-table:: Dimensional Stride Comparison
   :widths: 25 35 40
   :header-rows: 1

   * - Layout Model
     - X-Sweep Neighbor Stride
     - Y-Sweep Neighbor Stride
   * - **Row-Major**
     - $\pm 1$ (Coalesced)
     - $\pm \text{Stride}$ (Strided)
   * - **Column-Major**
     - $\pm \text{Height}$ (Strided)
     - $\pm 1$ (Coalesced)
     - **Tiled (32x32)**
     - Local: $\pm 1$ / Edge: Complex
     - Local: $\pm 32$ / Edge: Complex


L1/L2 Cache Hit Rate
=====================

Measured for each of the 9 (grid size, layout) conditions using the
compute sweep kernels.

Metrics
-------

.. list-table::
   :header-rows: 1

   * - Concept
     - Metric
     - Meaning
   * - L1 hit rate
     - ``l1tex__t_sector_hit_rate.pct``
     - Fraction of L1/texture cache accesses served from cache rather than
       going to L2/DRAM.
   * - L2 hit rate
     - ``lts__t_sector_hit_rate.pct``
     - Fraction of L2 cache accesses served from cache rather than going
       to DRAM.

Measurement
-----------

Profiled once per condition (access pattern, not timing, so no repeats
needed), then reduced to a per-(grid size, layout) average:

.. code-block:: bash

   ncu --nvtx --nvtx-include "Bench_500_RowMajor_0/" [...one --nvtx-include per condition...] \
       --metrics l1tex__t_sector_hit_rate.pct,lts__t_sector_hit_rate.pct \
       --export layout_bench_report --force-overwrite \
       ./build/benchmark_cuda RUNS --layout-benchmark

   ncu --import layout_bench_report.ncu-rep --csv --page raw > layout_bench.csv

.. code-block:: python

   # parse NVTX range -> GridDim/Layout, keep Sweep kernels, average per group
   summary = pd.read_csv("layout_bench.csv")[[..., 'l1tex__t_sector_hit_rate.pct',
                                               'lts__t_sector_hit_rate.pct']]
   summary[['GridDim', 'Layout']] = extract_from_nvtx_range(summary)  # Bench_<dim>_<layout>_<rep>
   agg = summary[summary.Kernel.str.contains('Sweep')] \
       .groupby(['GridDim', 'Layout'])[['L1_HitRate_pct', 'L2_HitRate_pct']].mean()
   agg.to_csv("hit_rate_avg_by_layout.csv")

Output: a 9-row table (``hit_rate_avg_by_layout.csv``) of average L1/L2
hit rate per (grid size, layout) condition.

Memory Coalescing and Achieved Bandwidth
=========================================

In addition to L1/L2 hit rates, two further metrics were collected for
each of the 9 (grid size, layout) conditions.

Metrics
-------

.. list-table::
   :header-rows: 1

   * - Concept
     - Metric
     - Meaning
   * - Load coalescing efficiency
     - ``smsp__sass_average_data_bytes_per_sector_mem_global_op_ld.pct``
     - Fraction of bytes fetched per sector that were actually useful for
       global loads. Low values mean threads in a warp access
       non-contiguous addresses, wasting sector bandwidth.
   * - Store coalescing efficiency
     - ``smsp__sass_average_data_bytes_per_sector_mem_global_op_st.pct``
     - Same as above, for global stores.
   * - Achieved DRAM bandwidth
     - ``dram__bytes.sum.per_second``
     - Bytes/second actually moved to/from DRAM during the kernel.
   * - DRAM throughput vs. peak
     - ``dram__throughput.avg.pct_of_peak_sustained_elapsed``
     - Achieved DRAM throughput as a percentage of the GPU's theoretical
       peak sustained bandwidth. Already a ratio; no further normalization
       needed.

Measurement
-----------

Same ``ncu``/NVTX setup as above, extended with four more metrics
(check availability first with ``ncu --query-metrics | grep -i coalesc``
/ ``dram__`` if any are rejected on your GPU):

.. code-block:: bash

   ncu --nvtx --nvtx-include "Bench_500_RowMajor_0/" [...] \
       --metrics l1tex__t_sector_hit_rate.pct,lts__t_sector_hit_rate.pct,\
   smsp__sass_average_data_bytes_per_sector_mem_global_op_ld.pct,\
   smsp__sass_average_data_bytes_per_sector_mem_global_op_st.pct,\
   dram__bytes.sum.per_second,dram__throughput.avg.pct_of_peak_sustained_elapsed \
       --export layout_bench_report_full --force-overwrite \
       ./build/benchmark_cuda RUNS --layout-benchmark

   ncu --import layout_bench_report_full.ncu-rep --csv --page raw > layout_bench_full.csv

.. code-block:: python

   # same parse/aggregate as before, plus convert DRAM bytes/sec -> GB/s
   summary['DRAM_BW_GBps'] = summary['DRAM_BW_Bps'] / 1e9
   agg = summary[summary.Kernel.str.contains('Sweep')] \
       .groupby(['GridDim', 'Layout'])[[
           'L1_HitRate_pct', 'L2_HitRate_pct', 'LoadCoalescing_pct',
           'StoreCoalescing_pct', 'DRAM_BW_GBps', 'DRAM_Throughput_pct']].mean()
   agg.to_csv("full_metrics_by_layout.csv")

Output: a 9-row table (``full_metrics_by_layout.csv``) combining cache hit
rate, coalescing efficiency, and achieved bandwidth per condition.

Runtime Performance Benchmark (Wall-Clock Timing)
==================================================

The sections above characterize *why* a layout performs the way it does,
using hardware counters captured through ``ncu``. This section instead
captures *how fast* each of the 9 (grid size, layout) conditions actually
runs, using the benchmark binary's own timers with no profiler attached.

Metrics
-------

.. list-table::
   :header-rows: 1

   * - Metric
     - Meaning
   * - ``AvgConversionSec``
     - Average time spent reformatting data into the target layout (e.g.
       re-indexing into column-major or tiled order) before the compute
       sweeps run. Zero-cost for the baseline layout the data already
       lives in, non-trivial for Tiled 32x32.
   * - ``AvgComputeSec``
     - Average time spent executing the finite-volume sweep kernels
       themselves, excluding layout conversion.
   * - ``AvgTotalSec``
     - Average end-to-end time per run (conversion + compute + any
       remaining fixed overhead).
   * - ``PaddingPct``
     - Percentage of the allocated grid occupied by padding/ghost-cell
       overhead rather than active domain data. Depends on tile size and
       grid dimensions, so it can differ between layouts at the same
       grid size.

Measurement
-----------

No profiler attached — the binary times itself, averaged over 10 repeats
per condition to reduce timing noise:

.. code-block:: bash

   ./build/benchmark_cuda 10 --layout-benchmark | tee layout_timing_sweep.txt

.. code-block:: text

   GridX  GridY  Layout       AvgConversionSec  AvgComputeSec  AvgTotalSec  PaddingPct
   500    500    RowMajor     0.000309767       3.91442e-05   0.000348911  0.8016
   ...

.. code-block:: python

   # regex-parse the console table into columns, cast numeric, save
   timing = parse_fixed_width("layout_timing_sweep.txt",
       columns=["GridX", "GridY", "Layout", "AvgConversionSec",
                "AvgComputeSec", "AvgTotalSec", "PaddingPct"])
   timing.to_csv("layout_timing_sweep.csv", index=False)

Output: a 9-row table (``layout_timing_sweep.csv``) of average conversion
time, compute time, total time, and padding overhead per condition —
the data plotted in ``images/Speed.png`` below.

Cache/Bandwidth Metrics vs. Runtime Timing
-------------------------------------------

.. list-table:: What Each Benchmark Captures
   :widths: 25 35 40
   :header-rows: 1

   * - 
     - Cache / Coalescing / Bandwidth
     - Runtime Timing Sweep
   * - Tool
     - ``ncu`` (Nsight Compute profiler)
     - Direct wall-clock timers in the benchmark binary
   * - Granularity
     - Per-kernel hardware counters
     - Per-run aggregate (conversion, compute, total)
   * - Columns
     - ``L1_HitRate_pct``, ``L2_HitRate_pct``, ``LoadCoalescing_pct``,
       ``StoreCoalescing_pct``, ``DRAM_BW_GBps``, ``DRAM_Throughput_pct``
     - ``AvgConversionSec``, ``AvgComputeSec``, ``AvgTotalSec``,
       ``PaddingPct``
   * - Question answered
     - *Why* is a layout fast or slow (cache/DRAM behavior)?
     - *How fast* is it in absolute wall-clock time, including any
       layout-conversion overhead?
   * - Repeats
     - 1 (access pattern is deterministic, not timing-sensitive)
     - 10 (averaged to reduce timing noise)

Because Tiled 32x32 requires an explicit re-layout pass, its
``AvgConversionSec`` is roughly 20-50x higher than Row-/Column-Major at
every grid size tested, and its ``PaddingPct`` is also higher whenever the
grid dimensions are not exact multiples of the tile size. Row-Major and
Column-Major track each other closely on conversion and total time, with
the larger differences showing up in ``AvgComputeSec`` instead, consistent
with the coalescing behavior described above.

.. image:: images/Speed.png
   :alt: Different computing averages for the memory Layouts
   :width: 600px
   :height: 400px
   :scale: 50%
   :align: center

.. image:: images/L1_vs._L2.png
   :alt: L1 vs. L2 cache hit rate and achieved bandwidth for the memory Layouts
   :width: 600px
   :height: 400px
   :scale: 50%
   :align: center
