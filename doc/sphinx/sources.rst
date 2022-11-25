
Sources and transfer terms
==========================

External sources are

-  :math:`S_n =` Source of plasma ions

-  :math:`S_p =` Source of pressure, related to energy source
   :math:`S_E = \frac{3}{2}S_p`

In the simulations carried out so far, these source functions are both
constant between midplane and X-point, and zero from X-point to target.

Transfer channels
-----------------

There are several transfer channels and sinks for particles, energy and
momentum due to rates of recombination, ionisation, charge exchange,
electron-neutral excitation, and elastic collisions with units of
m\ :math:`^{-3}`\ s\ :math:`^{-1}`:

.. math::

   \begin{aligned}
     \mathcal{R}_{rc} &=& n^2\left<\sigma v\right>_{rc}   \qquad \mbox{\textrm{(Recombination)}} \\
     \mathcal{R}_{iz} &=&  nn_n\left<\sigma v\right>_{iz} \qquad \mbox{\textrm{(Ionisation)}} \\
     \mathcal{R}_{cx} &=& nn_n\left<\sigma v\right>_{cx} \qquad \mbox{\textrm{(Charge exchange)}} \\
     \mathcal{R}_{el} &=& nn_n\left<\sigma v\right>_{el} \qquad \mbox{\textrm{(Elastic collisions)}}\end{aligned}

where :math:`n` is the plasma density; :math:`n_n` is the neutral gas
density; :math:`\sigma_{cx}` is the cross-section for charge exchange;
:math:`\sigma_{rc}` is the cross-section for recombination; and
:math:`\sigma_{iz}` is the cross-section for ionisation. Each of these
processes’ cross-section depends on the local density and temperatures,
and so changes in time and space as the simulation evolves.

-  :math:`S =` Net recombination i.e neutral source (plasma particle
   sink). Calculated as Recombination - Ionisation:

   .. math::

      \begin{aligned}
        S &=& \mathcal{R}_{rc} - \mathcal{R}_{iz}\end{aligned}

-  :math:`R =` Cooling of the plasma due to radiation, and plasma
   heating due to 3-body recombination at temperatures less than 5.25eV.

   .. math::

      \begin{aligned}
          R &=& \left(1.09 T_e - 13.6\textrm{eV}\right)\mathcal{R}_{rc} \qquad \mbox{\textrm{(Recombination)}}\\
          &+& E_{iz}\mathcal{R}_{iz}  \qquad \mbox{\textrm{(Ionisation)}} \\
          &+& \left(1\textrm{eV}\right)\mathcal{R}_{ex} \qquad \mbox{\textrm{(Excitation)}} \\
          &+& R_{z,imp} \qquad \mbox{\textrm{(Impurity radiation)}}
        \end{aligned}

   The factor of 1.09 in the recombination term, together with factor of
   :math:`3/2` in :math:`E` below, is so that recombination becomes a
   net heat source for the plasma at :math:`13.6 / 2.59 = 5.25`\ eV.
   :math:`E_{iz}` is the average energy required to ionise an atom,
   including energy lost through excitation.

   If excitation is not included (``excitation = false``) then following
   Togo *et al.*, :math:`E_{iz}` is chosen to be 30eV. If excitation is
   included, then :math:`E_{iz}` should be set to :math:`13.6`\ eV.

-  :math:`E =` Transfer of energy to neutrals.

   .. math::

      \begin{aligned}
          E &=& \frac{3}{2} T_e \mathcal{R}_{rc} \qquad \mbox{\textrm{(Recombination)}} \\
          &-& \frac{3}{2} T_n \mathcal{R}_{iz}  \qquad \mbox{\textrm{(Ionisation)}} \\
          &+& \frac{3}{2}\left(T_e - T_n\right)\mathcal{R}_{cx} \qquad \mbox{\textrm{(Charge exchange)**}} \\
          &+& \frac{3}{2}\left(T_e - T_n\right)\mathcal{R}_{el} \qquad \mbox{\textrm{(Elastic collisions)**}}
        \end{aligned}

   (**) Note that if the neutral temperature is not evolved, then
   :math:`T_n = T_e` is used to calculate the diffusion coefficient
   :math:`D_n`. In that case, :math:`T_n` is set to zero here, otherwise
   it would cancel and leave no CX energy loss term.

-  :math:`F =` Friction, a loss of momentum from the ions, due to charge
   exchange and recombination. The momentum of the neutrals is not
   currently modelled, so instead any momentum lost from the ions is
   assumed to be transmitted to the walls of the machine.

   .. math::

      \begin{aligned}
        F &=& m_iV_{||}\mathcal{R}_{rc} \qquad \mbox{\textrm{(Recombination)}} \\
        &-& m_iV_{||n}\mathcal{R}_{iz} \qquad \mbox{\textrm{(Ionisation)}} \\
        &+& m_i\left(V_{||} - V_{||n}\right)\mathcal{R}_{cx}  \qquad \mbox{\textrm{(Charge exchange)}} \\
        &+& m_i\left(V_{||} - V_{||n}\right)\mathcal{R}_{el}  \qquad \mbox{\textrm{(Elastic collisions)}}\end{aligned}

All transfer channels are integrated over the cell volume using
Simpson’s rule:

.. math:: S = \frac{1}{6J_C}\left( J_LS_L + 4J_CS_C + J_RS_R \right)

where :math:`J` is the Jacobian of the coordinate system, corresponding
to the cross-section area of the flux tube, and subscripts :math:`L`,
:math:`C` and :math:`R` refer to values at the left, centre and right of
the cell respectively.

Recycling
---------

The flux of ions (and neutrals) to the target plate is recycled and
re-injected into the simulation. The fraction of the flux which is
re-injected is controlled by ``frecycle``:

::

   frecycle = 0.95     # Recycling fraction

The remaining particle flux (5% in this case) is assumed to be lost from
the system. Note that if there are any external particle sources, then
this fraction must be less than 1, or the number of particles in the
simulation will never reach steady state.

Of the flux which is recycled, a fraction ``fredistribute`` is
redistributed along the length of the domain, whilst the remainder is
recycled at the target plate

::

   fredistribute = 0.8  # Fraction of recycled neutrals redistributed evenly along length

The weighting which determines how this is redistributed is set using
``redist_weight``:

::

   redist_weight = h(y - pi)  # Weighting for redistribution

which is normalised in the code so that the integral is always 1. In
these expressions :math:`y` is uniform in cell index, going from
:math:`0` to :math:`2\pi` between the boundaries. The above example
therefore redistributes the neutrals evenly (in cell index) from
half-way along the domain to the end.

When neutrals are injected, some assumptions are needed about their
energy and momentum

-  When redistributed, neutrals are assumed to arrive with no net
   parallel momentum (so nothing is added to :math:`NV_n`), and they are
   assumed to have the Franck-Condon energy (3.5eV currently)

-  When recycled from the target plate, neutrals are assumed to have a
   parallel momentum away from the target, with a thermal speed
   corresponding to the Franck-Condon energy, and is also added to the
   pressure equation. NOTE: This maybe should be one or the other, but
   not both...
