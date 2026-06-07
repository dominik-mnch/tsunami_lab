# Running Tsunami Lab on Draco

This guide describes how to log in to the Draco HPC cluster, copy or clone this
project, download the Chapter 6 input data, compile the solver on Draco, and run
interactive and batch jobs with Slurm.

Replace placeholders such as `<USER>` and `<REPO_URL>` with your own URZ user
name and repository URL.

## 1. Log in to Draco

Draco is reachable through the login nodes. If you are outside the university
network, connect to the university VPN first.

```bash
ssh <USER>@login1.draco.uni-jena.de
# or
ssh <USER>@login2.draco.uni-jena.de
```

The login nodes are shared machines. Use them for file transfers, editing,
compiling, and submitting jobs. Do not run long or expensive simulations on a
login node; use Slurm allocations for that.

Useful first checks:

```bash
hostname
pwd
sinfo -s
module avail
```

## 2. Choose a working directory

Use `$HOME` for source code and small files. For larger simulation input and
output, prefer a work file system such as `/work` if it is available to your
account.

```bash
mkdir -p ~/tsunami_runs
cd ~/tsunami_runs

# Optional: use /work for larger runs if your account has access.
mkdir -p /work/$USER/tsunami_runs
```

## 3. Get the source code onto Draco

### Option A: clone the Git repository on Draco

This is usually the cleanest option.

```bash
cd ~/tsunami_runs
git clone <REPO_URL> tsunami_lab
cd tsunami_lab
```

If the Catch2 submodule is not included after cloning, initialize submodules:

```bash
git submodule update --init --recursive
```

### Option B: copy the local workspace with scp

From a local PowerShell terminal on Windows:

```powershell
scp -r "G:\VSCodeProjekte\tsunami_lab" <USER>@login1.draco.uni-jena.de:~/tsunami_runs/
```

Then log in and enter the copied directory:

```bash
ssh <USER>@login1.draco.uni-jena.de
cd ~/tsunami_runs/tsunami_lab
```

For repeated uploads, `sftp`, WinSCP, FileZilla, or a Git workflow is more
comfortable than copying the whole directory every time.

## 4. Download the Chapter 6 input data

The assignment provides the input archive through the university cloud. Download
it directly on Draco so the large data does not have to pass through your local
machine.

```bash
cd ~/tsunami_runs/tsunami_lab
mkdir -p data/ch6_input

wget https://cloud.uni-jena.de/s/CqrDBqiMyKComPc/download/data_in.tar.xz \
  -O tsunami_lab_data_in.tar.xz

# Inspect the archive layout first.
tar -tf tsunami_lab_data_in.tar.xz | head

# Extract into a dedicated input directory.
tar -xf tsunami_lab_data_in.tar.xz -C data/ch6_input

# Check which netCDF and point-of-interest files are available.
find data/ch6_input -maxdepth 3 -type f | sort | head -50
```

Use the paths shown by `find` as arguments to `tsunami2d`. If the archive
contains an additional top-level directory, keep that directory in the path.

If you already have the data locally and want to upload it instead, use:

```powershell
scp -r "G:\VSCodeProjekte\tsunami_lab\data\output" <USER>@login1.draco.uni-jena.de:~/tsunami_runs/tsunami_lab/data/
```

## 5. Load software modules

This project is built with SCons, C++17, and netCDF. SCons is not available by
default on Draco, so load Python and install it into your user account. On Draco,
`nc-config` may already be available without loading a netCDF module; if it is,
you do not need an extra netCDF module.

```bash
module purge

# Python is needed for SCons.
module load tools/python/3.8
pip install --user scons

# Load a C++ compiler if needed.
module avail gcc

# Example compiler module; adapt it to the names shown by module avail.
module load compiler/gcc/11.3.0

# Make sure user-installed Python tools are on PATH.
export PATH="$HOME/.local/bin:$PATH"

which scons
which g++
which nc-config
nc-config --version
```

The build uses `nc-config --libs --cflags`, so `nc-config` must be available
before running SCons. If `which nc-config` fails, search for a netCDF module and
load it as a fallback:

```bash
module avail 2>&1 | grep -i netcdf
module spider netcdf
module load <NETCDF_MODULE_NAME>
which nc-config
```

## 6. Compile the solver on Draco

Build from the repository root:

```bash
cd ~/tsunami_runs/tsunami_lab
scons -j 8 mode=release
```

Expected executables:

```bash
ls -lh build/tsunami_lab build/tests
```

Optional sanity check:

```bash
./build/tests
```

If the build fails with `nc-config: command not found`, load a netCDF module and
run SCons again. If `scons` is not found, check `$HOME/.local/bin` and the
Python module.

## 7. Run a short interactive test

Request an interactive Slurm allocation on the `short` partition:

```bash
salloc --partition=short --ntasks=1 --cpus-per-task=1 --time=00:15:00
```

Draco forwards the shell to the allocated compute node after the allocation is
granted. Confirm that you are no longer on a login node:

```bash
hostname
```

Run a small test case:

```bash
cd ~/tsunami_runs/tsunami_lab
./build/tsunami_lab 100 1 0 10000 0 1 1 1 2.5 1d outflow dam_break_1d 10 0 5 0 5000
```

For a small two-dimensional tsunami test with project-local example data:

```bash
./build/tsunami_lab 1000m -200000 2500000 -750000 750000 1 1 60 \
  2d outflow tsunami2d \
  data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
  data/output/tohoku_gebco20_ucsb3_250m_displ.nc
```

For Chapter 6 data, replace the two netCDF paths with the matching files under
`data/ch6_input`.

Leave the compute node when finished:

```bash
exit
```

## 8. Submit a batch job

Interactive jobs are useful for quick checks. Production runs should use batch
jobs. Create `run_tsunami_draco.sbatch` in the repository root:

```bash
nano run_tsunami_draco.sbatch
```

Paste this script and adapt the input file paths and runtime:

```bash
#!/bin/bash
#SBATCH --job-name=tsunami_lab
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --hint=nomultithread
#SBATCH --time=03:00:00
#SBATCH --output=logs/tsunami_lab.%j.out
#SBATCH --error=logs/tsunami_lab.%j.err

set -euo pipefail

module purge
module load tools/python/3.8
module load compiler/gcc/11.3.0
export PATH="$HOME/.local/bin:$PATH"

which nc-config
nc-config --version

cd "$HOME/tsunami_runs/tsunami_lab"
mkdir -p logs solutions stations/output

echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"

srun ./build/tsunami_lab 1000m -200000 2500000 -750000 750000 1 1 7200 \
  2d outflow tsunami2d \
  data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
  data/output/tohoku_gebco20_ucsb3_250m_displ.nc

echo "Job finished at: $(date)"
```

If `nc-config` is not available by default in the batch job environment, add the
same netCDF fallback module that made `which nc-config` work during compilation.

Submit and monitor the job:

```bash
mkdir -p logs
sbatch run_tsunami_draco.sbatch
squeue --me
```

After completion, inspect the output files. Replace `<JOBID>` with the ID shown
by `sbatch`.

```bash
cat logs/tsunami_lab.<JOBID>.out
cat logs/tsunami_lab.<JOBID>.err
sacct -j <JOBID> -X
```

Cancel a pending or running job if needed:

```bash
scancel <JOBID>
```

## 9. Reproduce the Chapter 6 scenarios

Use the same command-line arguments and input files that were used for the
Chapter 6 local runs. The executable syntax is:

```bash
./build/tsunami_lab NX NY X_LOWER X_UPPER Y_LOWER Y_UPPER K SOLVER_MODE END_TIME PROPAGATION [PROP_PARAMS] SETUP [SETUP_PARAMS]
./build/tsunami_lab <RES>m X_LOWER X_UPPER Y_LOWER Y_UPPER K SOLVER_MODE END_TIME PROPAGATION [PROP_PARAMS] SETUP [SETUP_PARAMS]
```

Typical tsunami setup:

```bash
./build/tsunami_lab 250m X_LOWER X_UPPER Y_LOWER Y_UPPER K 1 END_TIME \
  2d outflow tsunami2d \
  path/to/bathymetry.nc \
  path/to/displacement.nc
```

Keep one directory per scenario so outputs do not overwrite each other:

```bash
mkdir -p runs/tohoku_250m
cd runs/tohoku_250m
../../build/tsunami_lab 250m ... 2d outflow tsunami2d ...
```

Compare final output files, station output, and plots against the Chapter 6
local results.

## 10. Measure cluster performance

The assignment asks for time per cell and iteration while excluding setup and
file I/O. In this project, the relevant numerical work happens in the
`while( l_simTime < l_endTime )` loop in `src/main.cpp`, but that loop also
writes station CSV data and periodic netCDF output. For the assignment metric,
the solver now prints `time stepping seconds`, `time steps`, and `time per cell
and iteration` at the end of each run.

For each run, record:

- `nx`: number of cells in x-direction
- `ny`: number of cells in y-direction
- `iterations`: number of time steps
- `time_stepping_seconds`: wall-clock time spent only in the stepping loop

Then compute:

```text
time_per_cell_iteration = time_stepping_seconds / (nx * ny * iterations)
```

Run the same scenario locally and on Draco with the same resolution, solver,
end time, and compiler optimization mode. Report both values and the speedup:

```text
speedup = local_time_per_cell_iteration / draco_time_per_cell_iteration
```

Use batch jobs for final measurements and keep the `sbatch` script, JobID,
stdout, stderr, and exact command-line arguments with the results.

## 11. Useful Draco commands

```bash
# Cluster and partition overview
sinfo
sinfo -s
sinfo -o "%20N %10c %10m %35f %10G"

# Job handling
sbatch run_tsunami_draco.sbatch
squeue --me
sacct -j <JOBID> -X
scancel <JOBID>

# Interactive allocation
salloc --partition=short --ntasks=1 --cpus-per-task=1 --time=00:15:00

# Modules
module avail
module list
module purge
module load <MODULE_NAME>

# Data transfer from local machine
scp -r <LOCAL_PATH> <USER>@login1.draco.uni-jena.de:<REMOTE_PATH>
sftp <USER>@login1.draco.uni-jena.de
```

## 12. Common problems

`scons: command not found`

Load Python and install SCons for your user:

```bash
module load tools/python/3.8
pip install --user scons
export PATH="$HOME/.local/bin:$PATH"
```

`nc-config: command not found`

First check whether it is available without a module. If not, search the full
module tree and load the module that provides it:

```bash
which nc-config
module avail 2>&1 | grep -i netcdf
module spider netcdf
module load <NETCDF_MODULE_NAME>
which nc-config
```

`undefined reference to std::filesystem...`

With GCC 12.2, `std::filesystem` is normally part of the standard C++ library.
Draco's system library stack may still require the compatibility library
`stdc++fs` explicitly at link time, especially when `nc-config` adds system
library paths such as `-L/usr/lib64`. The build script links `stdc++fs` after
the netCDF libraries. Pull the latest build script, clean stale objects, and
rebuild:

```bash
git pull
scons -c
scons -j 8 mode=release
```

Also verify that SCons sees the expected compiler:

```bash
which gcc
gcc --version
which g++
g++ --version
scons --debug=presub -j 1 mode=release
```

If the error remains with GCC 12.2 after a clean rebuild, inspect the final link
command printed by `scons --debug=presub`; it should use `g++`, not `gcc`, and
the link command should contain `-lstdc++fs` after `-lnetcdf`:

```bash
scons --debug=presub -j 1 mode=release
```

If needed, force the compiler explicitly:

```bash
scons -c
scons -j 8 mode=release CXX=g++
```

Job stays pending

Check the reason with `squeue --me -l`. Shorter runtimes, fewer resources, or a
different partition may schedule faster.

Simulation accidentally runs on login node

Stop it and request a Slurm allocation with `salloc` or submit through `sbatch`.
Use `hostname` to check where the shell is running.