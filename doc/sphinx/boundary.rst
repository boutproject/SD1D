
Boundary conditions
===================

Upstream: Symmetry
------------------

Symmetry boundary conditions are applied at the upstream side,
corresponding to zero flow through the boundary.

.. math:: \partial_{||}n = 0 \qquad \partial_{||}p = 0 \qquad \partial_{||}T_e = 0 \qquad V_{||} = 0 \qquad nV_{||} = 0

Since the boundary is half-way between grid points, this is implemented
as

.. math::

   \begin{aligned}
   n_0 &=& n_1 \\
   p_0 &=& p_1 \\
   T_{e,0} &=& T_{e,1} \\
   V_{||,0} &=& -V_{||,1} \\
   nV_{||,0} &=& -nV_{||,1}\end{aligned}

Downstream: Sheath
------------------

Boundary conditions are applied to the velocity and the heat flux:

-  At the left boundary a no-flow condition is applied:

   .. math::

      \begin{aligned}
          V_{||} &=& 0 \nonumber \\
          \partial_{||}T_e &=& 0 \nonumber
        \end{aligned}

-  At the right boundary is a sheath boundary:

   .. math::

      \begin{aligned}
          V_{||} &\ge& v_s \nonumber \\
          \partial_{||}T_e &=& 0 \nonumber
        \end{aligned}

   where the inequality is implemented by switching from a Dirichlet to
   a Neumann boundary if :math:`V_{||} > v_s` in front of the boundary.

   The critical speed into the sheath, :math:`v_s` is sensitive to
   assumptions on the thermodynamics of the sheath, taking the
   form: [1]_

   .. math:: v_s = \left( \frac{e\left(T_e + \gamma T_i\right)}{m_i}\right)^{1/2}

   where :math:`T_e` is the electron temperature (in eV), :math:`T_i` is
   the ion temperature, :math:`\gamma` is the ratio of specific heats.
   For isothermal flow :math:`\gamma=1`, for adiabatic flow with
   isotropic pressure :math:`\gamma=5/3`, and for one-dimensional
   adiabatic flow :math:`\gamma=3`. Here we are assuming
   :math:`T_e = T_i` and :math:`\partial_{||}T_e` so take the isothermal
   case. This therefore becomes:

   .. math::

      v_s = \left( \frac{p}{n}\right)^{1/2}
          \label{eq:sheath_speed}

Note: If the sheath velocity is subsonic, then waves can propagate in
from the boundary. Their domain of dependence is outside the simulation
domain, so these waves can cause numerical instabilities.

Several boundary conditions are available for the density and pressure,
including free boundaries and Neumann (zero gradient). These are
controlled by settings ``density_sheath`` and ``pressure_sheath``.
Density can have the following values:

0. Free boundary, linearly extrapolating the value from inside the
   domain

   .. math:: n_{-1} = 2n_{-2} - n_{-3}

1. Neumann (zero gradient)

   .. math:: n_{-1} = n_{-2}

2. Constant flux

   .. math:: n_{-1/2} = n_{-2}v_{-2}J_{-2} / \left( v_s J_{-1/2} \right)

   where the Jacobian factors :math:`J` account for a changing flux tube
   cross-section area.

Pressure can have the following values:

0. Free boundary, linearly extrapolating the value from inside the
   domain

   .. math:: p_{-1} = 2p_{-2} - p_{-3}

1. Neumann (zero gradient)

   .. math:: p_{-1} = p_{-2}

2. Constant energy flux :math:`\frac{5}{2}pv + \frac{1}{2}nv^3`

   .. math:: 5 p_{-1/2} = \left( 5p_{-2}v_{-2} + n_{-2}v_{-2}^3\right) / v_s - n_{-1/2}v_s^2
