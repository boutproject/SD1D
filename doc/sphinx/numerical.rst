
Numerical methods
=================

All variables are defined at the same location (collocated). Several
different numerical methods are implemented, to allow testing of their
accuracy and robustness.

Advection terms :math:`\nabla\cdot\left(\mathbf{b}V_{||}f\right)`
-----------------------------------------------------------------

.. _`sec:fluxsplit`:

Flux splitting, MinMod limiter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default method uses a combination of HLL-style flux splitting and
MinMod slope limiting. Terms of the form
:math:`\nabla\cdot\left(\mathbf{b} f\right)` are implemented as fluxes
through cell boundaries:

.. math:: \nabla\cdot\left(\mathbf{b} V f\right)_i \simeq \frac{1}{J\Delta y} \left[ F_{i+1/2} - F_{i-1/2}\right]

where :math:`F` is the flux. This is calculated by linearly
interpolating the velocity to the cell edges

.. math:: V_{i+1/2} = \frac{1}{2}\left(V_{i} + V_{i+1}\right)

The field being advected, :math:`f`, is reconstructed from the cell
centre values :math:`f_i` onto cell edges :math:`f^L_i` and
:math:`f^R_i`:

.. math:: f^L_i = f_i - \frac{1}{2}s \qquad f^R_i = f_i + \frac{1}{2}s

where the slope :math:`s` is limited using the MinMod method:

.. math::

   s = \left\{\begin{array}{ll}
   0 & \textrm{if $\operatorname{sign}(f_{i+1} - f_{i}) \neq \operatorname{sign}(f_{i} - f_{i-1})$} \\
   f_{i+1} - f_{i} & \textrm{if $\left|f_{i+1} - f_{i}\right| < \left|f_{i} - f_{i-1}\right|$} \\
   f_{i} - f_{i-1} & \textrm{otherwise}
   \end{array}\right.

In order to handle waves travelling both left and right, flux splitting
handles characteristics moving left differently from characteristics
moving right. In general this is problem dependent and computationally
expensive, so here we adopt a simple approximation similar to an HLL
splitting [2]_. We assume that the fastest waves in the system travel
with speed :math:`a` (the sound speed) with respect to the flow, so that
there are waves travelling with :math:`V+a` and :math:`V-a`. If the flow
speed is supersonic then these waves are only in one direction, but for
subsonic flows there is a flux in both directions. The fluxes between
cells are calculated using:

.. math::

   F_{i+1/2} = \left\{\begin{array}{ll}
   f^R_iV_{i+1/2} & \textrm{if $V_{i+1/2} > a$} \\
   f^L_{i+1}V_{i+1/2} & \textrm{if $V_{i+1/2} < -a$} \\
   f^R_i\frac{1}{2}\left(V_{i+1/2} +a\right) + f^L_{i+1}\frac{1}{2}\left(V_{i+1/2} - a\right) & \textrm{otherwise}
   \end{array}\right.
   \label{eq:splitfluxes}

Hence for subsonic flows the flux becomes
:math:`V_{i+1/2}\frac{1}{2}\left(f^R_i + f^L_{i+1}\right) + \frac{a}{2}\left(f^R_i - f^L_{i+1}\right)`,
where the second term is a diffusion. When the solution is smooth,
:math:`f^R_{i}\simeq f^L_{i+1}`, the numerical method becomes central
differencing and the diffusion goes to zero as :math:`\Delta x^2`.
Oscillatory solutions introduce dissipation, and the method becomes
increasingly upwind as the flow becomes sonic.

.. _`sec:nonlinflux`:

Nonlinear fluxes
~~~~~~~~~~~~~~~~

When advecting quantities which are a nonlinear combination of
variables, such as :math:`nV_{||}`, conservation properties can be
slightly improved by using the following interpolation [3]_  [4]_  [5]_:

.. math:: \left(fg\right)^R = \frac{1}{2}\left(f^Rg^C + f^Cg^R\right)

where superscript :math:`C` denotes cell centre, and :math:`R` right
hand side. This method is implemented, using MinMod interpolation for
each variable.

.. _`sec:central`:

Central differencing
~~~~~~~~~~~~~~~~~~~~

Central difference schemes have an advantage over upwind schemes, in
that they do not need to take account of wave speeds. The simple central
differencing scheme produces large unphysical oscillations, due to the
decoupling of odd and even points in collocated schemes, but can
(usually) be stabilised by adding dissipation. It is implemented here
for comparison with other schemes.

.. _`sec:skewform`:

Skew symmetric central differencing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simple modification to the central differencing scheme improves
numerical stability, coupling nearby points [6]_  [7]_ The idea is to
split the divergence terms into a “skew-symmetric” form

.. math:: \nabla\cdot\left(\mathbf{b} V_{||} f\right) = \frac{1}{2}\left[ \nabla\cdot\left(\mathbf{b} V_{||} f\right) + V_{||}\mathbf{b}\cdot\nabla f + f\nabla\cdot\left(\mathbf{b} V_{||}\right)\right]

Each of the terms on the right are then discretised with standard
2nd-order central differences. This method can avoid the need for
additional dissipation, or be stabilised with a smaller viscosity than
the simple central differencing method.

.. _`sec:viscos`:

Artificial viscosity
--------------------

Artificial viscosity (``viscos`` input) is implemented as a diffusion of
momentum in index space, so that the diffusion coefficient varies as
:math:`\Delta y^2`.

.. math:: \frac{\partial}{\partial t}\left(nV_{||}\right)_i = \ldots + \nu\left[ \left(V_{i+1} - V_i\right)J_{i+1/2} - \left(V_i - V_{i-1}\right)J_{i-1/2} \right]/J_i

where :math:`J` is the Jacobian, subscript :math:`i` indicates cell
index, and :math:`J_{i+1/2} = \left(J_i + J_{i+1}\right)/2`.

.. [1]
   K-U Riemann, J.Phys. D:Appl. Phys 24 (1991) 493-518

.. [2]
   A. Harten, P. D. Lax, and B. van Leer,”On Upstream Differencing and
   Godunov-Type Schemes for Hyperbolic Conservation Laws”, SIAM Review,
   25(1), pp. 35-61, 1983

.. [3]
   F.N.Felten, T.S.Lund “Kinetic energy conservation issues associated
   with the collocated mesh scheme for incompressible flow” J.Comp.Phys.
   215 (2006) 465-484

.. [4]
   F.N.Felten, T.S.Lund “Critical comparison of the collocated and
   staggered grid arrangements for incompressible turbulent flow” Report
   ADP013663

.. [5]
   Y.Morinishi et al. “Fully Conservative Higher Order Finite Difference
   Schemes for Incompressible Flow” J.Comp.Phys. 143 (1998) 90-124

.. [6]
   S.Pirozzoli “Stabilized non-dissipative approximations of Euler
   equations in generalized curvilinear coordinates” J.Comp.Phys. 230
   (2011) 2997-3014

.. [7]
   A.E.Honein, P.Moin “Higher entropy conservation and numerical
   stability of compressible turbulence simulations” J.Comp.Phys. 201
   (2004) 532-545
