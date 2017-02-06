
SD1D
====

SOL and Divertor in 1D. Simulates a plasma fluid in one dimension (along the magnetic field), interacting with a neutral gas fluid. 

* Evolves the density, momentum and pressure (internal energy) of both plasma and neutrals. 
* Includes exchange of particles, momentum and energy through ionisation, recombination and charge exchange
* Several different numerical methods are implemented, both upwind and central differencing

Author: Ben Dudson, University of York <benjamin.dudson@york.ac.uk>

Released under the GPL license

Installing BOUT++
-----------------

This version works with BOUT++ v4.0 or later

    git clone https://github.com/boutproject/BOUT-dev.git
    cd BOUT-dev

To run this model, preconditioning is strongly recommended, and requires the CVODE solver, part of [SUNDIALS](http://computation.llnl.gov/projects/sundials).
Tested with version 2.6.0. To enable CVODE, BOUT++ should be configured using

    ./configure --with-cvode

or

    ./configure --with-sundials

(which then also enables the IDA solver). Compile BOUT++ with

    make

Installing SD1D
---------------

Once BOUT++ is configured and compiled,  checkout a copy of SD1D

    git clone https://github.com/boutproject/SD1D.git

then compile SD1D, specifying the path to BOUT-dev

    cd SD1D
    make BOUT_TOP=/path/to/BOUT-dev/

Running examples
----------------

Once compiled, launch an MPI job to run the examples

    mpirun -np 4 ./sd1d -d case-01

where "-d case-01" specifies the directory to use for inputs and outputs.


