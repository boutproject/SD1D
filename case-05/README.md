Case 5: 90% recycling, density feedback
=======================================

Uses a PI controller to set the upstream density, with
a fixed power source and 90% recycling at the target

Note: This example sets the solver type to be "cvode",
which requires SUNDIALS. The "beuler" solver will also
work well with this, but the BOUT++ default "pvode" solver
will perform poorly because it lacks a good preconditioner.