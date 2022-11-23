Plasma model
============

Equations for the plasma density :math:`n`, pressure :math:`p` and
momentum :math:`m_inV_{||i}` are evolved:

.. math::

   \begin{aligned}
     \frac{\partial n}{\partial t} &=& - \nabla\cdot\left( \mathbf{b}V_{||} n\right) + S_n - S\\
     \frac{\partial}{\partial t}\left(\frac{3}{2}p\right) &=& -\nabla\cdot\mathbf{q} + V_{||}\partial_{||}p + S_p - E - R \\
     \frac{\partial}{\partial t}\left(m_i nV_{||}\right) &=& -\nabla\cdot\left(m_inV_{||}\mathbf{b}V_{||}\right) - \partial_{||} p - F\\
     j_{||} &=& 0 \\
     T_i &=& T_e = \frac{1}{2}\frac{p}{en} \\
     \mathbf{q} &=& \frac{5}{2}p\mathbf{b}V_{||} - \kappa_{||e}\partial_{||}T_e\end{aligned}

Which has a conserved energy:

.. math:: \int_V \left[ \frac{1}{2}m_inV_{||i}^2 + \frac{3}{2}p \right] dV

The heat conduction coefficient :math:`\kappa_{||e}` is a nonlinear
function of temperature :math:`T_e`:

.. math:: \kappa_{||e} = \kappa_0 T_e^{5/2}

where :math:`\kappa_0` is a constant. See
section `8 <#sec:heatconduction>`__ for details.

Operators are:

.. math:: \partial_{||}f = \mathbf{b}\cdot\nabla f \qquad \nabla_{||} f = \nabla\cdot\left(\mathbf{b} f\right)



.. _`sec:heatconduction`:

Heat conduction
---------------

Spitzer heat conduction is used

.. math:: \kappa_{||e} = 3.2\frac{ne^2T\tau_e}{m_e} \simeq 3.1\times 10^4 \frac{T^{5/2}}{\ln \Lambda}

which has units of W/m/eV so that in the formula
:math:`q = -\kappa_{||e} \nabla T`, :math:`q` has units of Watts per
m\ :math:`^2` and :math:`T` has units of :math:`eV`. This uses the
electron collision time:

.. math:: \tau_e = \frac{6\sqrt{2}\pi^{3/2}\epsilon_0^2\sqrt{m_e}T_e^{3/2}}{\ln \Lambda e^{2.5} n} \simeq 3.44\times 10^{11} \frac{T_e^{3/2}}{\ln \Lambda n}

in seconds, where :math:`Te` is in eV, and :math:`n` is in
m\ :math:`^{-3}`.

Normalising by the quantities in table `1 <#tab:normalisation>`__ gives

.. math:: \hat{\kappa}_{||e} = 3.2 \hat{n}\hat{T}_e\frac{m_i}{m_e}\tau_e\Omega_{ci}

where hats indicate normalised (dimensionless) variables.
