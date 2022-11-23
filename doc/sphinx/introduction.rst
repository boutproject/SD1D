This is a model of a 1D fluid, assuming equal ion and electron
temperatures, no electric fields or currents.

Getting started
===============

First get a copy of development branch of BOUT++. You can download a
tarball from https://github.com/boutproject/BOUT-dev, but it is strongly
recommended you use Git:

::

   $ git clone https://github.com/boutproject/BOUT-dev.git

Configure and make BOUT-dev, including SUNDIALS. This is available from
http://computation.llnl.gov/projects/sundials, and is needed for
preconditioning to work correctly.

::

   $ cd BOUT-dev
   $ ./configure --with-sundials
   $ make

The user manual for BOUT++ is in subdirectory of BOUT-dev called
"manual", and contains more detailed instructions on configuring and
compiling BOUT++. This will build the core library code, which is then
used in each model or test case (see the examples/ subdirectory)

Next download a copy of SD1D into the BOUT-dev/examples subdirectory.
This isn’t strictly necessary, but it makes the "make" command simpler
(otherwise you add an argument ``BOUT_TOP=/path/to/BOUT-dev/`` to make)

::

   BOUT-dev/examples/$ git clone https://github.com/boutproject/SD1D.git
   BOUT-dev/examples/$ cd SD1D
   BOUT-dev/examples/SD1D $ make

Hopefully you should see something like:

::

     Compiling  sd1d.cxx
     Compiling  div_ops.cxx
     Compiling  loadmetric.cxx
     Compiling  radiation.cxx
     Linking sd1d

Here the main code is in "sd1d.cxx" which defines a class with two
methods: ``init()``, which is run once at the start of the simulation to
initialise everything, and ``rhs()`` which is called every timestep. The
function of rhs() is to calculate the time derivative of each evolving
variable: In the ``init()`` function the evolving variables are added to
the time integration solver (around line 192). This time integration
sets the variables to a value, and then runs ``rhs()``. Starting line
782 of ``sd1d.cxx`` you’ll see the density equation, calculating
``ddt(Ne)``. ``Ne`` is the evolving variable, and ``ddt()`` is a
function which returns a reference to a variable which holds the
time-derivative of the given field.

BOUT++ contains many differential operators (see
``BOUT-dev/include/difops.hxx``), but work has been done on improving
the flux conserving Finite Volume implementations, and they’re not yet
in the public repository. These are defined in ``div_ops.hxx`` and
``div_ops.cxx``.

The atomic rates are used in ``sd1d.cxx`` starting around line 641, and
are defined in ``radiation.cxx`` and ``radiation.hxx``.

To run a simulation, enter:

::

   $ ./sd1d -d case-01

This will use the "case-01" subdirectory for input and output. All the
options for the simulation are in ``case-01/BOUT.inp``.

The output should be a whole bunch of diagnostics, printing all options
used (which also goes into log file BOUT.log.0), followed by the timing
for each output timestep:

::

   Sim Time  |  RHS evals  | Wall Time |  Calc    Inv   Comm    I/O   SOLVER

   0.000e+00          1       1.97e-02     1.9    0.0    0.2   21.6   76.3
   5.000e+03        525       1.91e-01    89.0    0.0    0.6    1.1    9.3
   1.000e+04        358       1.30e-01    88.8    0.0    0.6    1.4    9.2
   1.500e+04        463       1.68e-01    89.2    0.0    0.6    1.3    8.9
   2.000e+04        561       2.02e-01    89.6    0.0    0.6    1.1    8.7
   2.500e+04        455       1.65e-01    89.2    0.0    0.6    1.2    9.1

The simulation time (first column) is normalised to the ion cyclotron
frequency (as SD1D started life as part of a turbulence model), which is
stored in the output as "``Omega_ci``". So each output step is 5000 /
``Omega_ci`` = 104.4 microseconds. The number of internal timesteps is
determined by the solver, and determines the number of times the
``rhs()`` function was called, which is given in the second column. If
this number starts steadily increasing, it’s often a sign of numerical
problems.

To analyse the simulation, the data is stored in the "case-01"
subdirectory along with the input. You can use IDL or Python to look at
the "Ne", "NVi", "P" variables etc. which have the same names as in the
``sd1d.cxx`` code. See section `6 <#sec:output>`__ for details of the
output variables and their normalisation. The evolving variables should
each be 4D, but all dimensions are of size 1 except for the time and
parallel index (200). Please see the BOUT++ user manual for details of
setting up the Python and IDL reading ("collect") routines.

Examples
--------

Case 1: Without heat conduction (Euler’s equations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Removing heat conduction reduces the system to fluid (Euler) equations
in 1D. Note that in this case the boundary condition
(equation `[eq:sheath_speed] <#eq:sheath_speed>`__) is subsonic, because
the adiabatic fluid sound speed is

.. math:: c_s = \left( \frac{\gamma p}{n}\right)^{1/2} \qquad \gamma = 5/3

In this case the sources of particles and energy are uniform across the
grid.

Case 2: Localised source region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same equations are solved, but here the sources are only in the
first half of the domain, applied with a Heaviside function so the
sources abruptly change.

Case 3: Heat conduction
~~~~~~~~~~~~~~~~~~~~~~~

We now add Spitzer heat conduction, the :math:`\kappa_{||e}` term in the
pressure equation. This coefficient depends strongly on temperature, and
severely limits the timestep unless preconditioning is used. Here we use
the CVODE solver with preconditioning of the electron heat flux. In
addition to improving the speed of convergence, this preconditioning
also improves the numerical stability.

Case 4: Recycling, neutral gas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The plasma equations are now coupled to a similar set of equations for
the neutral gas density, pressure, and parallel momentum. A fixed
particle and power source is used here, and a 20% recycling fraction.
Exchange of particles, momentum and energy between neutrals and plasma
occurs through ionisation, recombination and charge exchange.

Case 5: High recycling, upstream density controller
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses a PI feedback controller to set the upstream density
to :math:`1\times 10^{19}`\ m\ :math:`^{-3}`. This adjusts the input
particle source to achieve the desired density, so generally needs some
tuning to minimise transient oscillations. This is controlled by the
inputs

::

   density_upstream = 1e19
   density_controller_p = 1e-2
   density_controller_i = 1e-3

The input power flux is fixed, specified in the input as
:math:`20`\ MW/m\ :math:`^2`:

::

   [P]    # Plasma pressure P = 2 * Ne * T
   powerflux = 2e7  # Input power flux in W/m^2

The recycling is set to 95%

::

   frecycle = 0.95

**NOTE**: This example is under-resolved; a realistic simulation would
use a higher resolution, but would take longer. To increase resolution
adjust ``ny``:

::

   ny = 200    # Resolution along field-line

Rather than 200, a more realistic value is about 600 or higher
with a uniform mesh. An alternative is to compress grid cells
closer to the target by varying the grid spacing ``dy``.

Non-uniform mesh
----------------

An example of using a non-uniform grid is in ``diffusion_pn``. The
location :math:`l` along the field line as a function of normalised cell
index :math:`y`, which goes from :math:`0` at the upstream boundary to
:math:`2\pi` at the target, is

.. math::

   l = L\left[ \left(2 - \delta y_{min}\right)\frac{y}{2\pi} -\left(1-\delta y_{min}\right)\left(\frac{y}{2\pi}\right)^2\right]
   \label{eq:nonuniform_l}

where :math:`0<\delta y_{min}<1` is a parameter which sets the size of
the smallest grid cell, as a fraction of the average grid cell size. The
grid cell spacing :math:`\delta y` therefore varies as

.. math:: \delta y = \frac{L}{N_y} \left[ 1 + \left(1-\delta y_{min}\right)\left(1-\frac{y}{\pi}\right)\right]

This is set in the BOUT.inp settings file, under the ``mesh`` section:

::

   dy = (length / ny) * (1 + (1-dymin)*(1-y/pi))

In order to specify the size of the source region, the normalised cell
index :math:`y` at which the location :math:`l` is a given fraction of
the domain length must be calculated. This is done by solving for
:math:`y` in equation `[eq:nonuniform_l] <#eq:nonuniform_l>`__.

.. math:: y_{xpt} = \pi\left[2 - \delta y_{min} - \sqrt{\left(2-\delta y_{min}\right)^2 - 4\left(1-\delta y_{min}\right) f_{source}}\right]/\left(1-\delta y_{min}\right)

which is calculated in the BOUT.inp file as

::

   y_xpt = pi * ( 2 - dymin - sqrt( (2-dymin)^2 - 4*(1-dymin)*source ) ) / (1 - dymin)

where ``source`` is the fraction :math:`f_{source}` of the length over
which the source is spread. This is then used to calculate sources,
given a total flux. For density:

::

   source = (flux/(mesh:source*mesh:length))*h(mesh:y_xpt - y)

which switches on the source for :math:`y < y_{xpt}` using a Heaviside
function, then divides the flux by the length of the source region
:math:`f_{source}L` to get the volumetric sources.

