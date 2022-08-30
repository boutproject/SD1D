
Neutral model
=============

The number of equations solved is controlled by the following parameters
in the input file:

::

   [NVn]
   evolve = true # Evolve neutral momentum?

   [Pn]
   evolve = true # Evolve neutral pressure? Otherwise Tn = Te model

Neutral density is always evolved, so turning off evolution of momentum
and pressure (setting both of the above to false) reduces the neutral
model to a simple diffusion model (next section). By turning on the
momentum equation

Diffusive model
---------------

In the simplest neutral model, neutral gas is modelled as a fluid with a
density :math:`n_n` which diffuses with a diffusion coefficient
:math:`D_n`:

.. math:: \frac{\partial n_n}{\partial t} = \nabla\cdot\left(D_n\nabla n_n\right) + S - n_n / \tau_n

The temperature of the neutrals is assumed to be the same as the ions
:math:`T_n = T_i`.Diffusion of neutrals depends on the neutral gas
temperature, and on the collision rate:

.. math:: D_n = v^2_{th,n} / \left(\nu_{cx} + \nu_{nn}\right)

where :math:`v_{th,n} = \sqrt{eT_n/m_i}` is the thermal velocity of a
neutral atom; :math:`\nu_{cx} = n\sigma_{cx}` is the charge-exchange
frequency, and :math:`\sigma_{nn} = v_{th,n} n_n a_0` is the
neutral-neutral collision frequency where
:math:`a_0 \simeq \pi \left(5.29\times 10^{-11}\right)^2` m\ :math:`^2`
is the cross-sectional area of a neutral Hydrogen atom. In order to
prevent divide-by-zero problems at low densities, which would cause
:math:`D` to become extremely large, the mean free path of the neutrals
is limited to :math:`1`\ m.

An additional loss term is required in order to prevent the particle
inventory of the simulations becoming unbounded in detached simulations,
where recycling no longer removes particles from the system. This
represents the residence time for neutral particles in the divertor
region, which in [Togo 2013] was set to around :math:`10^{-4}`\ s.

Neutral fluid model
-------------------

A more sophisticated neutrals model can be used, which evolves the
neutral gas momentum and energy:

.. math::

   \begin{aligned}
     \frac{\partial n_n}{\partial t} &=& - \nabla\cdot\left( \mathbf{b}V_{||n} n_n\right) + {\nabla\cdot\left(D_n\nabla n_n\right)} + S - n_n / \tau_n\\
     \frac{\partial}{\partial t}\left(\frac{3}{2}p_n\right) &=& -V_{||n}\partial_{||}p_n + \nabla\cdot\left(\kappa_n\nabla T_n\right) + \nabla\cdot\left(D_nT_n\nabla n_n\right) + E \\
     \frac{\partial}{\partial t}\left(m_i nV_{||n}\right) &=& -\nabla\cdot\left(m_inV_{||n}\mathbf{b}V_{||n}\right) - \partial_{||} p + F\\\end{aligned}

where :math:`\kappa_n` is the neutral gas heat conduction coefficient.
This is assumed to be

.. math:: \kappa_n = n_n v_{th,n}^2 / \left(\nu_{cx} + \nu_{nn}\right)

i.e. similar to :math:`D_n` for the diffusive neutral model, but with a
factor of :math:`n_n`.

Note that if the diffusion term :math:`D_n` is retained in the neutral
density (:math:`n_n`) equation, then a corresponding term is needed in
the pressure (:math:`p_n`) equation. To remove these terms, set
``dneut`` to zero in the input options, which will set :math:`D_n = 0`.

The density diffusion term should not be included if the momentum is
evolved, and so is switched off if this is the case. The continuity
equation for :math:`n_n` is exact once the flow is known, so the
diffusive flux should be contained in the flow velocity :math:`V_{||n}`.
To see where this comes from, assume an isothermal neutral gas:

.. math::

   \begin{aligned}
     \frac{\partial n_n}{\partial t} &=& - \nabla\cdot\left( \mathbf{b}V_{||n} n_n\right) + S - n_n / \tau_n\\
     \frac{\partial}{\partial t}\left(m_i nV_{||n}\right) &=& -\nabla\cdot\left(m_inV_{||n}\mathbf{b}V_{||n}\right) - eT_n\partial_{||} n_n + F\end{aligned}

Dropping the inertial terms reduces the momentum equation to

.. math:: eT_n\partial_{||} n_n = F = \nu m_i n_n \left(V_{||i} - V_{||n}\right)

where :math:`\nu` is a collision frequency of the neutrals with the
ions, due to charge exchange, recombination and ionisation (i.e.
:math:`\nu_{cx} + \nu_{nn}` as used in the calculation of diffusion
coefficient :math:`D_n`). This gives an equation for the neutral flow
velocity:

.. math:: V_{||n} = V_{||i} - \frac{eT_n}{m_i n_n \nu}\partial_{||} n_n = \frac{1}{n_n}\frac{v_{th,n}^2}{\nu}\partial_{||} n_n

where :math:`v_{th} = \sqrt{eT_n / m_i}` is the neutral thermal speed,
as used in the calculation of :math:`D_n`. This gives a flux of neutrals

.. math:: n_nV_{||n} = n_nV_{||i} - D_n\partial_{||} n_n

Hence the diffusive flux is included in the balance between pressure
gradients and friction in the momentum equation.
