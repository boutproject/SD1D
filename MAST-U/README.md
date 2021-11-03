MAST-Upgrade like input
=======================

Includes:

- A PI controller to set the upstream density, with
  a fixed power source and 99% recycling at the target

- Non-uniform grid, packed near the target

- The backward Euler solver, to find steady state solutions quickly
  Note: This requires PETSc to be linked to BOUT++ when compiling
  (CMake option BOUT_USE_PETSC)

- Nitrogen impurity (0.5%) using ADAS coronal model
