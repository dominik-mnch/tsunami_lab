Draft Individual Phase
======================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

The project will conclude with an individual phase. To be able to properly execute on a plan, we want to first
make a draft discussing what we want to do in the individual phase and some risks and problems that could arise.

Simulating boats and other vessels in a tsunami
-----------------------------------------------

When starting a tsunami simulation, we want to be able to include boats or other vessels in the simulation.
They could be stationary or moving through the domain. We want to be able to determine whether or not the boat
would be able to withstand the tsunami or if it would sink. To do this we would have to come up with a metric
to determine the damage to the boat. This could be based on the forces acting on the boat, the water height, or other factors.

To model the vessel itself, we can have it start a single point in the domain and use a regular function to describe its
movement across the domain. Then, for each time step, we can check whether the boat is in a cell that is affected by the tsunami.

A potential expansion to this could also be to have the boat react to the tsunami and try to move out of the way or, more
generally, if it is even possible to move out of the way. We could potentially calculate the optimal path
for the boat to take to get away from the tsunami.

Another expansion could be to have a quicker and smaller boat saving the crew of the bigger boat. This way the boat would be destroyed
but the crew would be saved. We would also have to calculate an optimal route for the smaller boat to take to get to the bigger boat and then get away from the tsunami,
and calculate, wether it is even possible to do so in time. 

Some problems that could arise are:

- How do we come up with the metric for the damage to the boat? This is a very arbitrary choice and it could be 
  difficult to make our metric reflect reality.

- How do we model the boat itself with its many properties? Boats with different shapes, masses, velocities and materials will 
  react differently to the tsunami and modelling that can be difficult.

- How difficult would it be to model a potential escape path for the boat? This is a problem that's completely different 
  from the rest of the project and could be challenging to implement.


Further parallelization using CUDA
----------------------------------

The last non-individual phase of the project will be to parallelize the code using OpenMP. This is a good start to
harness the power of parallelization. Since the different parts of the main loop to calculate the next time step are
largely independent of each other, this is a very natural thing to do. 

However, we can go even further and use CUDA to parellelize the code on a GPU and then run this code on a GPU cluster.
This would allow us to run larger simulations and get results faster. This is especially useful for tsunami 
simulations where we want to be able to run large simulations with a lot of data.

We already have some basic understanding of using CUDA from the course "Efficient Computing" where we got a basic 
introduction to CUDA and parallel programming on GPUs. 

Some problems that could arise are:

- Parallelizing the code using CUDA is not an easy task and using the library effectively can be difficult. Even if the code
  looks solid, one small error can cause the program to run very slowly or even crash. Debugging it can also be difficult
  since it doesn't run on our machines and we have to use the output from the cluster to do it.

- Currently, our time step calculates a net update for each edge of a cell. This means that it applies updates to multiple
  cells in a single iteration of the loop. If we were to parallelize this loop completely, we would have to make sure that no
  conflicting updates are applied to the same cell at the same time. We could use ``atomicAdd`` for this (but this could make the
  code slower) or we could change the structure of the code so every loop iteration is responsible for a single cell update.

