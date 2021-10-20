
SD1D
====

SOL and Divertor in 1D. Simulates a plasma fluid in one dimension (along the magnetic field), interacting with a neutral gas fluid. 

* Evolves the density, momentum and pressure (internal energy) of both plasma and neutrals. 
* Includes exchange of particles, momentum and energy through ionisation, recombination and charge exchange
* Several different numerical methods are implemented, both upwind and central differencing

Author: Ben Dudson, University of York <benjamin.dudson@york.ac.uk>

Released under the GPL license

Installing with CMake
---------------------

This is probably the most straightforward method to use now.  First
configure BOUT++ and SD1D. To use the default options and minimal
dependencies just run:

    $ cmake . -B build

Alternatively the CMake build can be customised: See the [BOUT++
documentation](https://bout-dev.readthedocs.io/en/latest/user_docs/installing.html#cmake)
for examples of using `cmake` arguments, or edit the compile options
interactively before building:

    $ ccmake . -B build

The best solver currently available is "beuler", which requires PETSc.
If you have PETSc installed, then BOUT++ can be built with it by
setting `BOUT_USE_PETSC` to `ON` i.e. on the cmake command line include
`-DBOUT_USE_PETSC=ON`. For more information on compiling BOUT++ with
PETSc see [the BOUT++
manual](https://bout-dev.readthedocs.io/en/latest/user_docs/advanced_install.html#petsc).
If PETSc is not available, then it is highly recommended to enable
SUNDIALS, by setting `BOUT_USE_SUNDIALS` to `ON`.

During configuration
[BOUT++](https://github.com/boutproject/BOUT-dev/) will be
automatically downloaded as a submodule, together with some
dependencies (NetCDF and FFTW are assumed to be installed already,
along with optional dependencies like SUNDIALS and PETSc if they are
requested).  Once configured, run build to compile BOUT++ and then
SD1D:

    $ cmake --build build


Installing with Autotools
-------------------------

To compile with the Makefile system, first install BOUT++.
This version of SD1D works with BOUT++ v4.3 or later:

    git clone -b v4.3.0-rc https://github.com/boutproject/BOUT-dev.git
    cd BOUT-dev

To run this model, preconditioning is strongly recommended, and
requires the CVODE solver, part of
[SUNDIALS](http://computation.llnl.gov/projects/sundials).  Tested
with version 2.6.0. To enable CVODE, BOUT++ should be configured using

    ./configure --with-cvode

or

    ./configure --with-sundials

(which then also enables the IDA solver). Compile BOUT++ with

    make

Once BOUT++ is configured and compiled, checkout a copy of SD1D

    git clone https://github.com/boutproject/SD1D.git

then compile SD1D, specifying the path to BOUT-dev

    cd SD1D
    make BOUT_TOP=/path/to/BOUT-dev/

Running examples
----------------

Once compiled, with either CMake or Autotools, launch an MPI job to run the examples

    mpirun -np 4 ./sd1d -d case-01

where "-d case-01" specifies the directory to use for inputs and outputs.


