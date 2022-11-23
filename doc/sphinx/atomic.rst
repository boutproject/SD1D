
Atomic cross sections
=====================

Cross sections are approximated with semi-analytic expressions, obtained
from E.Havlickova but of unknown origin. For the purposes of calculating
these cross-sections, any temperatures below 1eV are set to 1eV. The
charge exchange cross-section is approximated as:

.. math::

   \sigma_{iz} = \left\{\begin{array}{ll}
   10^{-14} T^{1/3} & \textrm{if $T \ge 1$eV} \\
   10^{-14} & \textrm{if $T < 1$eV} \end{array}\right.

with units of :math:`[\textrm{m}^3/\textrm{s}]`. Ionisation is
calculated as

.. math::

   \sigma_{cx} = \left\{\begin{array}{ll}
   5.875\times 10^{-12}\cdot T^{-0.5151} \cdot 10^{-2.563/\log_{10}T} & \textrm{if $T \ge 20$eV} \\
   10^{-6}\cdot T^{-3.054}\cdot 10^{-15.72\exp\left(-\log_{10}T\right) + 1.603\exp\left(-\log^2_{10}T\right)} & \textrm{if $1$eV $ < T < 20$eV} \\
   7.638\times 10^{-21} & \textrm{if $T \le 1$eV}\end{array}\right.

Recombination rates are calculated using a :math:`9\times 9` table of
coefficients so is not reproduced here.

.. figure:: hydrogen.pdf
   :alt: Cross-sections [Thanks to E.Havlickova and H.Willett]
   :name: fig:sigma

   Cross-sections [Thanks to E.Havlickova and H.Willett]

Plots of these cross-sections are shown in figure `1 <#fig:sigma>`__.
There are a few anomalies with this: charge exchange always has the
highest cross-section of any process, and ionisation has a jump at
:math:`20`\ eV. The ionisation and charge exchange rates do not depend
on density, but recombination does so a typical range of values is
shown.
