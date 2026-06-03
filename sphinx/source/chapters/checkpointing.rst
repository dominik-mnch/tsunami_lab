Checkpointing
=============

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Overview
--------

Long-running simulations — such as the Tohoku or Chile tsunami events at high resolution — may
take several hours or days to complete.  If the process is interrupted (power loss, job
scheduler timeout, user abort), all progress is lost without a recovery mechanism.

To address this, we added a **checkpointing** system that periodically saves the full simulation
state to disk.  When the simulator is started again with the same output directory, it detects
the checkpoint automatically and resumes from the last saved state without any additional
command-line arguments.


File Layout
-----------

The checkpointing system uses two files, both written to the ``solutions/`` directory:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - File
     - Purpose
   * - ``solutions/checkpoint.nc``
     - Simulation *parameters* (grid dimensions, cell sizes, origin, end time) and the most
       recently committed simulation time ``sim_time``.  Updated on every output step.
   * - ``solutions/solution.nc``
     - Time-series of the full wave field (``height``, ``momentum_x``, ``momentum_y``,
       ``bathymetry``).  A new frame is appended every 25 time steps.


Checkpoint File Format
----------------------

``checkpoint.nc`` stores only global attributes — no dimensions or variables:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Attribute
     - Type
     - Description
   * - ``nx``
     - ``int``
     - Number of cells in x-direction.
   * - ``ny``
     - ``int``
     - Number of cells in y-direction.
   * - ``dx``
     - ``float``
     - Cell width in x-direction (metres).
   * - ``dy``
     - ``float``
     - Cell width in y-direction (metres).
   * - ``origin_x``
     - ``float``
     - x-coordinate of the domain origin.
   * - ``origin_y``
     - ``float``
     - y-coordinate of the domain origin.
   * - ``end_time``
     - ``float``
     - Configured simulation end time (seconds).
   * - ``sim_time``
     - ``float``
     - Simulation time of the most recent written output frame (seconds).
       Updated atomically after every ``writeTimeStep`` call.


Write Strategy
--------------

Two levels of durability are maintained:

1. **Every output frame** (every 25 time steps): a full wave-field snapshot is appended to
   ``solution.nc`` via ``writeTimeStep``.  Immediately afterwards, ``nc_sync`` is called on
   both ``solution.nc`` and ``checkpoint.nc`` so the data survives a sudden process kill.

2. **Every time step**: only ``sim_time`` in ``checkpoint.nc`` is updated via
   ``overwriteCheckpointSimTime``.  This keeps the checkpoint time stamp as fresh as possible
   at minimal I/O cost, but since no corresponding wave-field frame is written, this value
   is only informational until the next full output.

The effective resume granularity is therefore one output interval (25 time steps).  If the
process is killed between two output frames, the checkpoint time in ``checkpoint.nc`` may be
ahead of the last frame in ``solution.nc``; the reader handles this correctly (see below).


Resume Detection
----------------

At startup, before parsing any command-line arguments, ``main`` checks for the existence of
``solutions/checkpoint.nc``:

.. code-block:: text

   if solutions/checkpoint.nc exists:
       read CheckpointData from checkpoint.nc + solution.nc
       if simTime > 0:
           set l_useCheckpoint = true
       else:
           fall through to normal (fresh) start

The condition ``simTime > 0`` guards against two edge cases:

* The checkpoint file was created but the simulation was killed before the first output frame
  was written (``solution.nc`` has zero time steps).  In this case ``simTime`` is set to ``0``
  and a fresh start is performed.
* A stale checkpoint from a fully completed run.  If ``sim_time == end_time``, the guard
  still passes, but the time loop ``while( simTime < endTime )`` exits immediately.

When ``l_useCheckpoint`` is true, the normal setup and propagation objects are replaced by a
``Checkpoint`` setup instance.  All grid parameters (``nx``, ``ny``, ``dx``, ``dy``,
``origin_x``, ``origin_y``, ``end_time``) are taken from the checkpoint file, overriding the
command-line values.  The existing ``solution.nc`` is kept and the ``NetCdf`` writer opens it
in append mode.


Checkpoint Setup
----------------

The ``Checkpoint`` class (``src/setups/Checkpoint/``) is a regular ``Setup`` implementation
that loads the last committed wave field from the checkpoint pair.  It is constructed with the
path to ``checkpoint.nc``; the sibling ``solution.nc`` is located automatically.

Internally it delegates all netCDF I/O to ``NetCdf::readCheckpoint``, which:

1. Opens ``checkpoint.nc`` and reads all global attributes.
2. Derives the path to ``solution.nc`` (same directory).
3. Queries the ``time`` dimension length in ``solution.nc``.
4. If zero frames exist, returns ``simTime = 0`` (no-resume sentinel).
5. Otherwise scans backwards for the last frame whose time equals ``end_time`` (the commit
   marker written at end of simulation), falling back to the last frame with a finite time
   value.
6. Reads ``height``, ``momentum_x``, ``momentum_y``, and ``bathymetry`` at that frame.

Cell lookup at query coordinates uses nearest-neighbour indexing:

.. math::

   i_x = \mathrm{clamp}\!\left(\left\lfloor \frac{x - x_0}{\Delta x} \right\rfloor,\, 0,\, n_x - 1\right)

   i_y = \mathrm{clamp}\!\left(\left\lfloor \frac{y - y_0}{\Delta y} \right\rfloor,\, 0,\, n_y - 1\right)

   \text{index} = i_y \cdot n_x + i_x


NetCdf Writer in Append Mode
----------------------------

When resuming, ``solution.nc`` already contains frames from the previous run.  The ``NetCdf``
constructor detects the existing file and, instead of creating a new one:

1. Opens it with ``NC_WRITE``.
2. Queries all variable IDs (``time``, ``height``, ``momentum_x``, ``momentum_y``,
   ``bathymetry``) via ``nc_inq_varid``.
3. Sets ``m_timeStep`` to the current length of the ``time`` dimension so subsequent
   ``writeTimeStep`` calls append after the existing frames.

A new ``checkpoint.nc`` is always created (``NC_CLOBBER``) so that its ``end_time`` reflects
the current run's configured end time.


Usage
-----

Checkpointing is **fully automatic** — no extra flags are needed.

**Starting a simulation normally:**

.. code-block:: bash

   ./build/tsunami_lab 250m -200000 2500000 -750000 750000 1 1 7200 2d outflow tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc

**Resuming after an interruption** (run the exact same command again):

.. code-block:: bash

   ./build/tsunami_lab 250m -200000 2500000 -750000 750000 1 1 7200 2d outflow tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc

The simulator will print:

.. code-block:: text

   checkpoint detected — resuming from t = <simTime>
     resuming simulation at t = <simTime>

and continue appending to the existing ``solutions/solution.nc``.

**Starting fresh** (discarding a previous run):

Delete the checkpoint and solution files before invoking the simulator:

.. code-block:: bash

   rm solutions/checkpoint.nc solutions/solution.nc
