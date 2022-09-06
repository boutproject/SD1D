
.. _`sec:output`:

Outputs
=======

Output quantities are normalised, with the normalisation factors stored
in the output files

.. container:: center

   .. container::
      :name: tab:normalisation

      .. table:: Normalisation quantities

         ============ =========== ================
         Name         Description Units
         ============ =========== ================
         ``Nnorm``    Density     m\ :math:`^{-3}`
         ``Tnorm``    Temperature eV
         ``Cs0``      Speed       m/s
         ``Omega_ci`` Time        1/s
         ``rho_s0``   Length      m
         ============ =========== ================

The following variables are stored in the output file if they are
evolved:

.. container:: center

   +---------+------------------+---------------------------------------+
   | Name    | Description      | Normalisation                         |
   +=========+==================+=======================================+
   | ``Ne``  | Plasma density   | ``Nnorm`` [:math:`m^{-3}`]            |
   +---------+------------------+---------------------------------------+
   | ``NVi`` | Plasma flux      | ``Nnorm``\ :math:`\times`\ ``Cs0``    |
   |         |                  | [:math:`m^{-2}s^{-1}`]                |
   +---------+------------------+---------------------------------------+
   | ``P``   | Plasma pressure  | e\ :math:`\times`\                    |
   |         |                  |  ``Nnorm``\ :math:`\times`\ ``Tnorm`` |
   +---------+------------------+---------------------------------------+
   | ``Nn``  | Neutral density  | ``Nnorm`` [:math:`m^{-3}`]            |
   +---------+------------------+---------------------------------------+
   | ``NVn`` | Neutral flux     | ``Nnorm``\ :math:`\times`\ ``Cs0``    |
   |         |                  | [:math:`m^{-2}s^{-1}`]                |
   +---------+------------------+---------------------------------------+
   | ``Pn``  | Neutral pressure | e\ :math:`\times`\                    |
   |         |                  |  ``Nnorm``\ :math:`\times`\ ``Tnorm`` |
   +---------+------------------+---------------------------------------+

The following rates and coefficients are also stored:

.. container:: center

   +----------------+-------------------------+-------------------------+
   | Name           | Description             | Normalisation           |
   +================+=========================+=========================+
   | ``S``          | Sink of plasma density  | ``Nnorm`` :math        |
   |                |                         | :`\times`\ ``Omega_ci`` |
   |                |                         | [m\ :math:`^{           |
   |                |                         | -3}`\ s\ :math:`^{-1}`] |
   +----------------+-------------------------+-------------------------+
   | ``F``          | Sink of plasma momentum | :math:`m_i\times`       |
   |                |                         | ``Nnorm``:math:`\times` |
   |                |                         | ``Cs0`` :math:`\times`  |
   |                |                         | ``Omega_ci``            |
   |                |                         | [Nm\ :math:`^{-3}`]     |
   +----------------+-------------------------+-------------------------+
   | ``R``          | Radiative loss of       | :math:`e\times`         |
   |                | energy                  | ``Nnorm``:math:`\times` |
   |                |                         | ``Tnorm``:math:`\times` |
   |                |                         | ``Omega_ci``            |
   |                |                         | [Wm\ :math:`^{-3}`]     |
   +----------------+-------------------------+-------------------------+
   | ``E``          | Sink of plasma energy   | :math:`e\times`         |
   |                |                         | ``Nnorm``:math:`\times` |
   |                |                         | ``Tnorm``:math:`\times` |
   |                |                         | ``Omega_ci``            |
   |                |                         | [Wm\ :math:`^{-3}`]     |
   +----------------+-------------------------+-------------------------+
   | ``kappa_epar`` | Plasma thermal          |                         |
   |                | conduction              |                         |
   +----------------+-------------------------+-------------------------+
   | ``Dn``         | Neutral diffusion       |                         |
   |                | coefficient             |                         |
   +----------------+-------------------------+-------------------------+
   | ``flux_ion``   | Flux of ions to target  |                         |
   +----------------+-------------------------+-------------------------+

Note that the ``R`` term is energy which is lost from the system, whilst
``E`` is energy which is transferred between plasma and neutrals. For
all transfer terms (``S``, ``F``, ``R``) a positive value means a
transfer from plasma to neutrals.

To diagnose atomic processes, turn on ``diagnose = true`` in the input
settings (this is the default). Additional outputs contain the
contributions from each atomic process. They have the same normalisation
factors as the corresponding (``S``, ``F``, ``R``) term.

.. container:: center

   ========= =====================================================
   Name      Description
   ========= =====================================================
   ``Srec``  Sink of plasma particles due to recombination
   ``Siz``   Sink of plasma particles due to ionisation (negative)
   ``Frec``  Sink of plasma momentum due to recombination
   ``Fiz``   Sink of plasma momentum due to ionisation
   ``Fcx``   Sink of plasma momentum due to charge exchange
   ``Fel``   Sink of plasma momentum due to elastic collisions
   ``Rrec``  Radiation loss due to recombination
   ``Riz``   Radiation loss due to ionisation (inc. excitation)
   ``Rzrad`` Radiation loss due to impurities
   ``Rex``   Radiation loss due to electron-neutral excitation
   ``Erec``  Sink of plasma energy due to recombination
   ``Eiz``   Sink of plasma energy due to ionisation
   ``Ecx``   Sink of plasma energy due to charge exchange
   ``Eel``   Sink of plasma energy due to elastic collisions
   ========= =====================================================
