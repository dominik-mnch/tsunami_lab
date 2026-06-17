.. Tsunami Lab documentation master file, created by
   sphinx-quickstart on Thu Apr  9 12:28:55 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Tsunami Lab documentation
=========================

This is the documentation for the Tsunami Lab project by Magdalena Schwarzkopf (207559) and Dominik Münch (208100).

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Overview <self>
   chapters/f_wave_solver
   chapters/finite_volume_discretization
   chapters/bathymetry_boundary_conditions
   chapters/two-dimensional_solver
   chapters/large_data_input_and_output
   chapters/tsunami_simulations
   chapters/checkpointing_and_coarse_output
   chapters/optimization
   chapters/parallelization
   chapters/draft_individual_phase
   chapters/plan_individual_phase
   chapters/plan_individual_phase_full

You can find the source code documentation (generated with Doxygen and hosted with GitHub Pages) `here <https://dominik-mnch.github.io/tsunami_lab/doxygen/>`_.

Build Process
=============

The project uses `SCons <https://scons.org/>`_ as its build tool and `Nix Shell <https://nixos.org/>`_ to provide a reproducible development environment.

To set this up yourself you can follow this quick tutorial.

**First time setup:**

Install nix-shell:

.. code-block:: bash

   curl --proto '=https' --tlsv1.2 -sSf -L https://install.determinate.systems/nix | sh -s -- install

After doing so, you have to restart your terminal.

**Local build process:**

Run the nix-shell (the necessary packages are automatically included using the ``shell.nix`` file:

.. code-block:: bash

   nix-shell

You can now use SCons to compile the code:

.. code-block:: bash

   scons
   
To execute tests or the program itself you can use:

.. code-block:: bash

   ./build/tests
   ./build/tsunami_lab [nx] [ny] [x_lower] [x_upper] [y_lower] [y_upper] [k] [solver_mode] [end_time] [propagation] [prop_params] [setup] [setup_params]

Continuous Integration:
=======================

There are two GitHub actions pipelines set up that are triggered every time a commit is pushed to the ``main`` branch.

* ``main.yml``: A pipeline that runs a static code analysis using cppcheck, unit tests, sanitizer builds and Valgrind memory checks along with release builds.
* ``docs.yml``: A pipeline that builds Doxygen Code Documentation and Sphinx Documentation and deploys the resulting html to GitHub Pages.

