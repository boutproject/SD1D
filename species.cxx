
#include "species.hxx"
#include "div_ops.hxx"

FluidSpecies::FluidSpecies(std::string name, Options *opt, Solver *solver,
                           bool restarting,
                           BoutReal Nnorm, BoutReal Tnorm,
                           BoutReal Omega_ci, BoutReal Cs0)
  : name(name), Nnorm(Nnorm), Tnorm(Tnorm), Omega_ci(Omega_ci), Cs0(Cs0) {
  
  OPTION(opt, gamma_sound, 5. / 3); // Ratio of specific heats
  OPTION(opt, bndry_flux_fix, false);

  // Sheath boundary
  OPTION(opt, density_sheath, 0);   // Free boundary
  OPTION(opt, pressure_sheath, 0);  // Free boundary
  OPTION(opt, sheath_gamma, 6.5);   // Sheath heat transmission

  // Plasma anomalous transport
  OPTION(opt, anomalous_D, -1);
  OPTION(opt, anomalous_chi, -1);
  
  OPTION(opt, viscos, -1);            // Parallel viscosity
  
  OPTION(opt, ion_viscosity, false);  // Braginskii parallel viscosity

  BoutReal Coulomb = 6.6 - 0.5 * log(Nnorm * 1e-20) + 1.5 * log(Tnorm);
  tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3. / 2));
  
  solver->add(N, ("N" + name).c_str());
  solver->add(P, ("P" + name).c_str());
  solver->add(NV, ("NV" + name).c_str());
}

void FluidSpecies::evolve(BoutReal UNUSED(time)) {
  TRACE("Species::evolve");

  // Communicate evolving variables
  mesh->communicate(N, NV, P);

  Coordinates *coord = mesh->coordinates();

  // Floor small values
  P = floor(P, 1e-10);
  N = floor(N, 1e-10);

  Field3D Nlim = floor(N, 1e-5);

  V = NV / N; // Velocity
  T = P / N;  // Temperature
  
  {
    TRACE("Upper Y boundary (sheath)");

    ddt(P) = 0.0; // Need to set heat flux

    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      int jz = 0;

      // Outward flow velocity to >= Cs
      BoutReal Vout = sqrt(2.0 * T(r.ind, mesh->yend, jz));

      if (V(r.ind, mesh->yend, jz) > Vout) {
        // If plasma is faster than sound speed, go to plasma velocity
        Vout = V(r.ind, mesh->yend, jz);
      }

      BoutReal Nout;
      switch (density_sheath) {
      case 0: {
        // Free boundary on density (constant gradient)
        Nout = 0.5 *
               (3. * N(r.ind, mesh->yend, jz) - N(r.ind, mesh->yend - 1, jz));
        break;
      }
      case 1: {
        // Zero gradient
        Nout = N(r.ind, mesh->yend, jz);
        break;
      }
      case 2: {
        // Zero gradient particle flux N*Vi* J*dx*dz
        // Since Vi increases into the sheath, density should drop
        Nout =
            N(r.ind, mesh->yend, jz) * coord->J(r.ind, mesh->yend) *
            V(r.ind, mesh->yend, jz) /
            (0.5 *
             (coord->J(r.ind, mesh->yend) + coord->J(r.ind, mesh->yend + 1)) *
             Vout);
        break;
      }
      default:
        throw BoutException("Unrecognised density_sheath option");
      }
      
      // Prevent Nout from going negative
      // -> Flux is always to the wall
      if (Nout < 0.0)
        Nout = 0.0; 

      // Flux of particles is Ne*Vout
      BoutReal flux = Nout * Vout;

      BoutReal Pout;

      switch (pressure_sheath) {
      case 0: {
        // Free boundary  (constant gradient)
        Pout = 0.5 *
               (3. * P(r.ind, mesh->yend, jz) - P(r.ind, mesh->yend - 1, jz));
        break;
      }
      case 1: {
        // Zero gradient
        Pout = P(r.ind, mesh->yend, jz);
        break;
      }
      case 2: {
        // Use energy flux conservation to set pressure
        // (5/2)Pv + (1/2)nv^3 = const
        //
        Pout =
            ((5. * P(r.ind, mesh->yend, jz) * V(r.ind, mesh->yend, jz) +
              N(r.ind, mesh->yend, jz) * pow(V(r.ind, mesh->yend, jz), 3)) /
                 Vout -
             Nout * Vout * Vout) /
            5.;
        break;
      }
      default:
        throw BoutException("Unrecognised pressure_sheath option");
      }

      if (Pout < 0.0)
        Pout = 0.0;

      
      // Additional cooling
      BoutReal q = (sheath_gamma - 6) * T(r.ind, mesh->yend, jz) * flux;
      
      // Multiply by cell area to get power
      BoutReal heatflux =
        q *
        (coord->J(r.ind, mesh->yend) + coord->J(r.ind, mesh->yend + 1)) /
        (sqrt(coord->g_22(r.ind, mesh->yend)) +
         sqrt(coord->g_22(r.ind, mesh->yend + 1)));
      
      // Divide by volume of cell, and 2/3 to get pressure
      ddt(P)(r.ind, mesh->yend, jz) -=
        (2. / 3) * heatflux /
        (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));
      
      
      // Set boundary half-way between cells
      for (int jy = mesh->yend + 1; jy < mesh->LocalNy; jy++) {
        
        // velocity
        V(r.ind, jy, jz) = 2. * Vout - V(r.ind, mesh->yend, jz);

        // density
        N(r.ind, jy, jz) = 2 * Nout - N(r.ind, mesh->yend, jz);

        // NV. This can be negative, so set this to the flux
        // going out of the domain (zero gradient)
        NV(r.ind, jy, jz) = Nout * Vout;
        
        // temperature zero gradient (Neumann)
        T(r.ind, jy, jz) = T(r.ind, mesh->yend, jz);

        P(r.ind, jy, jz) = 2. * Pout - P(r.ind, mesh->yend, jz);
      }
    }
  }

  {
    TRACE("Lower Y boundary (no flow)");
    // No-flow boundary condition on lower Y boundary
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        for (int jy = 0; jy < mesh->ystart; jy++) {
          T(r.ind, jy, jz) = T(r.ind, mesh->ystart, jz);
          N(r.ind, jy, jz) = N(r.ind, mesh->ystart, jz);
          P(r.ind, jy, jz) = P(r.ind, mesh->ystart, jz);
          V(r.ind, jy, jz) = -V(r.ind, mesh->ystart, jz);
          NV(r.ind, jy, jz) = -NV(r.ind, mesh->ystart, jz);
        }
      }
    }
  }

  Field3D a = sqrt(gamma_sound * T); // Local sound speed

  {
    TRACE("Density");

    ddt(N) = -Div_par_FV_FS(N, V, a, bndry_flux_fix);

    if (volume_source) {
      ddt(Ne) += NeSource; // External volume source
    }

    if (anomalous_D > 0.0) {
      ddt(Ne) += Div_par_diffusion(anomalous_D, Ne);
    }
  }
  {
    TRACE("Momentum");
    ddt(NV) = -Div_par_FV_FS(NV, V, a, bndry_flux_fix) // Momentum flow
              - Grad_par(P);

    if (viscos > 0.) {
      ddt(NVi) += viscos * Div_par_diffusion_index(Vi);
    }
    
    if (anomalous_D > 0.0) {
      ddt(NVi) += Div_par_diffusion(anomalous_D * Vi, Ne);
    }
    
    if (ion_viscosity) {
      // Braginskii ion viscosity
      Field3D tau_i = sqrt(2 * mi_me) * Omega_ci * tau_e0 * pow(T, 1.5) / N;
      eta_i = (4. / 3) * 0.96 * N * tau_i * T;
      eta_i.applyBoundary("neumann");

      ddt(NVi) += Div_par_diffusion(eta_i, Vi);
    }
  }
  {
    TRACE("Pressure");

    // Note: ddt(P) set earlier for sheath
    ddt(P) += -Div_par_FV_FS(P, V, a, bndry_flux_fix) // Advection
             - (2. / 3) * P * Div_par(V)             // Compression
      ;
    
    // External source of energy
    ddt(P) += PeSource;
    
    // if (heat_conduction) {
    //   ddt(P) += (2. / 3) * Div_par_diffusion_upwind(kappa_epar, Te);
    // }

    if (anomalous_D > 0.0) {
      ddt(P) += Div_par_diffusion(anomalous_D * T, N);
    }
    
    if (anomalous_chi > 0.0) {
      ddt(P) += Div_par_diffusion(anomalous_chi, T);
    }
    
  }

  // Switch off evolution at very low densities
  for (auto i : ddt(N).region(RGN_NOBNDRY)) {
    if ((N[i] < 1e-5) && (ddt(N)[i] < 0.0)) {
      ddt(N)[i] = 0.0;
      ddt(NV)[i] = 0.0;
      ddt(P)[i] = 0.0;
    }
  }
}
