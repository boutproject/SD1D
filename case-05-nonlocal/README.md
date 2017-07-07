Case 5 (nonlocal): 90% recycling, density feedback
==================================================

Uses a PI controller to set the upstream density, with
a fixed power source and 90% recycling at the target

This version uses non-local electron heat transport.
This is enabled by setting "nonlocal_conduction = true"
in the [SD1D] section of BOUT.inp. 

See settings in the [non_local_parallel] section for settings
needed to get this to work correctly.

Two new outputs are "qe_nonlocal" and "qe_local" which are
the heat flux in normalised units from the two different calculations.



