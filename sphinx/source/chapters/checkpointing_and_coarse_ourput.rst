Checkpointing and Coarse Output
==================================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Coarse Output 
-------------

The coarse output is a downsampled version of the full resolution output. It allows for reduced data storage and processing requirements 
without having to deal with the full resolution data, which can be very large. 

The downsampling factor is determined by the parameter `K`, which specifies how many grid cells in each direction are combined into one cell 
in the coarse output. For example, if `K=2`, then each cell in the coarse output represents a 2x2 block of cells in the full resolution output.
If `K=1`, the coarse output is identical to the full resolution output.

Height, Bathymetry, and Velocity are averaged over the corresponding blocks of cells in the full resolution output to compute the values for the coarse output.
The coarse output is written to a separate NetCDF file, which can be used for visualization or further analysis.

Implementation details
----------------------

The coarse output implementation requires two pieces of work in the application code.

In `src/main.cpp`, the coarse factor `K` must be parsed from the command line and passed to the NetCDF writer. 
The main program already must now read `K` as part of the general input arguments, and it also needs to ensure that 
the domain grid dimensions are compatible with the coarse factor (for example, `NX % K == 0` and `NY % K == 0`), 
and then forward `K` into the `NetCdf` constructor.

In `src/io/NetCdf.h` / `src/io/NetCdf.cpp`, the constructor takes an additional parameter `i_k` and stores it as `m_k`. 
That value is used to define the coarse output dimensions (`nx/k` and `ny/k`) and to average blocks of `k*k` cells for 
height, bathymetry, and momentum. 
When writing each snapshot, the coarse writer iterates over the full-resolution field in steps of `k`, 
averages each block, and stores the aggregated values into the reduced output arrays.
The coarse output is written to a separate NetCDF file, which can be used for visualization or further analysis.

Simulating Tohoku earthquake and tsunami with coarse output
------------------------------------------------------------


./build/tsunami_lab 50m -200000 2500000 -750000 750000 50 1 3600 2d outflow tsunami2d \
    data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
    data/output/tohoku_gebco20_ucsb3_250m_displ.nc


Indiviual Contributions
-----------------------