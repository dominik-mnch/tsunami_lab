============================================================
CUDA Memory Layout Transformations: Column-Major and Tiled
============================================================

Memory Layout Frameworks
========================

To maximize global memory throughput on the GPU, data layout must align with hardware 
coalescing behaviors. Below are the three memory access models evaluated:

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

Profiling command
------------------

One repeat per condition is profiled (cache hit rate is a property of the
access pattern, not of run-to-run timing noise):

.. code-block:: bash

   ncu --nvtx \
       --nvtx-include "Bench_500_RowMajor_0/" \
       --nvtx-include "Bench_500_ColumnMajor_0/" \
       --nvtx-include "Bench_500_Tile32_0/" \
       --nvtx-include "Bench_1000_RowMajor_0/" \
       --nvtx-include "Bench_1000_ColumnMajor_0/" \
       --nvtx-include "Bench_1000_Tile32_0/" \
       --nvtx-include "Bench_4000_RowMajor_0/" \
       --nvtx-include "Bench_4000_ColumnMajor_0/" \
       --nvtx-include "Bench_4000_Tile32_0/" \
       --metrics l1tex__t_sector_hit_rate.pct,lts__t_sector_hit_rate.pct \
       --export layout_bench_report \
       --force-overwrite \
       ./build/benchmark_cuda RUNS --layout-benchmark

   ncu --import layout_bench_report.ncu-rep --csv --page raw > layout_bench.csv

Analysis
--------

.. code-block:: python

   import pandas as pd

   df = pd.read_csv("layout_bench.csv")

   cols = [
       'thread Domain:Push/Pop_Range:PL_Type:PL_Value:CLR_Type:Color:Msg_Type:Msg',
       'Kernel Name',
       'l1tex__t_sector_hit_rate.pct',
       'lts__t_sector_hit_rate.pct',
   ]

   summary = df[cols].copy()
   summary.columns = ['NVTX_Range', 'Kernel', 'L1_HitRate_pct', 'L2_HitRate_pct']

   summary = summary.dropna(subset=['NVTX_Range'])
   summary['NVTX_Range'] = summary['NVTX_Range'].str.extract(r'(Bench_\d+_\w+?_\d+)')
   extracted = summary['NVTX_Range'].str.extract(r'Bench_(\d+)_(\w+?)_(\d+)')
   summary['GridDim'] = extracted[0]
   summary['Layout']  = extracted[1]

   summary['L1_HitRate_pct'] = pd.to_numeric(summary['L1_HitRate_pct'], errors='coerce')
   summary['L2_HitRate_pct'] = pd.to_numeric(summary['L2_HitRate_pct'], errors='coerce')

   compute_kernels = summary[summary['Kernel'].str.contains('Sweep', na=False)]

   agg = (compute_kernels
          .groupby(['GridDim', 'Layout'])[['L1_HitRate_pct', 'L2_HitRate_pct']]
          .mean()
          .reset_index())

   agg.to_csv("hit_rate_avg_by_layout.csv", index=False)
   summary.to_csv("hit_rate_per_kernel.csv", index=False)

The result is a 9-row table (``hit_rate_avg_by_layout.csv``) with average
L1/L2 hit rate per (grid size, layout) condition.

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

Profiling command
------------------

.. code-block:: bash

   ncu --nvtx \
       --nvtx-include "Bench_500_RowMajor_0/" \
       --nvtx-include "Bench_500_ColumnMajor_0/" \
       --nvtx-include "Bench_500_Tile32_0/" \
       --nvtx-include "Bench_1000_RowMajor_0/" \
       --nvtx-include "Bench_1000_ColumnMajor_0/" \
       --nvtx-include "Bench_1000_Tile32_0/" \
       --nvtx-include "Bench_4000_RowMajor_0/" \
       --nvtx-include "Bench_4000_ColumnMajor_0/" \
       --nvtx-include "Bench_4000_Tile32_0/" \
       --metrics l1tex__t_sector_hit_rate.pct,lts__t_sector_hit_rate.pct,\
   smsp__sass_average_data_bytes_per_sector_mem_global_op_ld.pct,\
   smsp__sass_average_data_bytes_per_sector_mem_global_op_st.pct,\
   dram__bytes.sum.per_second,\
   dram__throughput.avg.pct_of_peak_sustained_elapsed \
       --export layout_bench_report_full \
       --force-overwrite \
       ./build/benchmark_cuda RUNS --layout-benchmark

   ncu --import layout_bench_report_full.ncu-rep --csv --page raw > layout_bench_full.csv

If a metric name is rejected, check what is available on the GPU first:

.. code-block:: bash

   ncu --query-metrics | grep -i coalesc
   ncu --query-metrics | grep -i "dram__"

Analysis
--------

.. code-block:: python

   import pandas as pd

   df = pd.read_csv("layout_bench_full.csv")

   cols = [
       'thread Domain:Push/Pop_Range:PL_Type:PL_Value:CLR_Type:Color:Msg_Type:Msg',
       'Kernel Name',
       'l1tex__t_sector_hit_rate.pct',
       'lts__t_sector_hit_rate.pct',
       'smsp__sass_average_data_bytes_per_sector_mem_global_op_ld.pct',
       'smsp__sass_average_data_bytes_per_sector_mem_global_op_st.pct',
       'dram__bytes.sum.per_second',
       'dram__throughput.avg.pct_of_peak_sustained_elapsed',
   ]

   summary = df[cols].copy()
   summary.columns = ['NVTX_Range', 'Kernel', 'L1_HitRate_pct', 'L2_HitRate_pct',
                       'LoadCoalescing_pct', 'StoreCoalescing_pct',
                       'DRAM_BW_Bps', 'DRAM_Throughput_pct']

   summary = summary.dropna(subset=['NVTX_Range'])
   summary['NVTX_Range'] = summary['NVTX_Range'].str.extract(r'(Bench_\d+_\w+?_\d+)')
   extracted = summary['NVTX_Range'].str.extract(r'Bench_(\d+)_(\w+?)_(\d+)')
   summary['GridDim'] = extracted[0]
   summary['Layout']  = extracted[1]

   for col in ['L1_HitRate_pct', 'L2_HitRate_pct', 'LoadCoalescing_pct',
               'StoreCoalescing_pct', 'DRAM_BW_Bps', 'DRAM_Throughput_pct']:
       summary[col] = pd.to_numeric(summary[col], errors='coerce')

   summary['DRAM_BW_GBps'] = summary['DRAM_BW_Bps'] / 1e9

   compute_kernels = summary[summary['Kernel'].str.contains('Sweep', na=False)]

   agg = (compute_kernels
          .groupby(['GridDim', 'Layout'])[[
              'L1_HitRate_pct', 'L2_HitRate_pct',
              'LoadCoalescing_pct', 'StoreCoalescing_pct',
              'DRAM_BW_GBps', 'DRAM_Throughput_pct'
          ]]
          .mean()
          .reset_index())

   agg.to_csv("full_metrics_by_layout.csv", index=False)
   summary.to_csv("full_metrics_per_kernel.csv", index=False)

The result is a 9-row table (``full_metrics_by_layout.csv``) combining cache
hit rate, coalescing efficiency, and achieved bandwidth for every
(grid size, layout) condition.