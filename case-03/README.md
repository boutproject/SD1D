Case 3: Plasma fluid with heat conduction
=========================================

Includes parallel heat conduction, with preconditioning
to speed up timestepping.

Stongly recommended to use the CVODE solver (part of SUNDIALS)
rather than PVODE, as this allows preconditioning.
To do this, BOUT++ should be configured with CVODE. See BOUT++ manual
for details.

