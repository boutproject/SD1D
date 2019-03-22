
#include "species.hxx"
#include "bout/constants.hxx"
#include "div_ops.hxx"
#include "field_factory.hxx"

using bout::globals::mesh;
using bout::globals::dump;

FluidSpecies::FluidSpecies(std::string name, Options *opt, Solver *solver,
                           Datafile &restart, bool restarting, BoutReal Nnorm,
                           BoutReal Tnorm, BoutReal Omega_ci, BoutReal Cs0)
    : name(name), Nnorm(Nnorm), Tnorm(Tnorm), Omega_ci(Omega_ci), Cs0(Cs0) {

  OPTION(opt, AA, 2.0); // Ion mass
  OPTION(opt, ZZ, 1.0); // Ion charge

  BoutReal rho_s0 = Cs0 / Omega_ci;

  OPTION(opt, gamma_sound, 5. / 3); // Ratio of specific heats
  OPTION(opt, bndry_flux_fix, false);

  // Sheath boundary
  OPTION(opt, sheath_outflow, true); // Out-flowing sheath boundary
  OPTION(opt, density_sheath, 0);  // Free boundary
  OPTION(opt, pressure_sheath, 0); // Free boundary
  OPTION(opt, sheath_gamma, 6.5);  // Sheath heat transmission

  // Plasma anomalous transport
  OPTION(opt, anomalous_D, -1);   // Input in m^2/s
  OPTION(opt, anomalous_chi, -1); // Input in m^2/s

  // Anomalous transport
  if (anomalous_D > 0.0) {
    // Normalise
    anomalous_D /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous D_perp = %e\n", anomalous_D);
  }
  if (anomalous_chi > 0.0) {
    // Normalise
    anomalous_chi /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous chi_perp = %e\n", anomalous_chi);
  }

  OPTION(opt, viscos, -1);           // Parallel viscosity
  OPTION(opt, ion_viscosity, false); // Braginskii parallel viscosity

  BoutReal Coulomb = 6.6 - 0.5 * log(Nnorm * 1e-20) + 1.5 * log(Tnorm);
  tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3. / 2));

  solver->add(N, ("N" + name).c_str());
  solver->add(P, ("P" + name).c_str());
  solver->add(NV, ("NV" + name).c_str());
  
  // Volume sources of particles and energy

  std::string source_string;
  FieldFactory ffact(mesh);

  Options *optne = Options::getRoot()->getSection("N" + name);
  optne->get("source", source_string, "0.0");
  NeSource = ffact.create2D(source_string, optne);

  Options *optpe = Options::getRoot()->getSection("P" + name);
  optpe->get("source", source_string, "0.0");
  PeSource = ffact.create2D(source_string, optpe);
  dump.addOnce(PeSource, "P" + name + "source");

  PeSource *= 0.5; /// <- Note: This is temporary, applies if Te = Ti
  
  // Normalise sources
  NeSource /= Nnorm * Omega_ci;
  PeSource /= SI::qe * Nnorm * Tnorm * Omega_ci;

  // Density controller
  OPTION(opt, density_upstream, -1); // Fix upstream density? [m^-3]
  if (density_upstream > 0.0) {
    // Fixing density
    density_upstream /= Nnorm;

    // Controller
    OPTION(opt, density_controller_p, 1e-2);
    OPTION(opt, density_controller_i, 1e-3);
    OPTION(opt, density_integral_positive, false);
    OPTION(opt, density_source_positive, true);

    density_error_lasttime = -1.0; // Signal no value

    // Save and load error integral from file, since
    // this determines the source function
    restart.add(density_error_integral, "density_error_integral");

    if (!restarting) {
      density_error_integral = 0.0;

      // Set density_error_integral so that
      // the input source is used
      density_error_integral = 1. / density_controller_i;
    }
    dump.addRepeat(NeSource, "N" + name + "source");
    NeSource0 = NeSource; // Save initial value
  } else {
    dump.addOnce(NeSource, "N" + name + "source");
  }
}

void FluidSpecies::evolve(BoutReal time) {
  TRACE("Species::evolve");

  // Communicate evolving variables
  mesh->communicate(N, NV, P);

  // Note: Here we are not using yup/ydown fields
  // and just writing boundary conditions into the one field
  N.mergeYupYdown();
  NV.mergeYupYdown();
  P.mergeYupYdown();
  
  Coordinates *coord = mesh->getCoordinates();

  // Floor small values
  P = floor(P, 1e-10);
  N = floor(N, 1e-10);

  V = NV / (AA * N); // Velocity. Note atomic mass factor AA
  T = P / N;         // Temperature

  /////////////////////////////////////////////////////
  // Boundaries

  {
    TRACE("Upper Y boundary (sheath)");

    ddt(P) = 0.0; // Need to set heat flux

    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      int jz = 0;

      // Out-flow velocity 
      BoutReal Vout = 0.0;
      if (sheath_outflow) {
        // Outward flow velocity to >= Cs
        // Note: Need to account for the species mass
        Vout = sqrt(2.0 * T(r.ind, mesh->yend, jz) / AA);
        
        if (V(r.ind, mesh->yend, jz) > Vout) {
          // If plasma is faster than sound speed, go to plasma velocity
          Vout = V(r.ind, mesh->yend, jz);
        }
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
      if (Nout < 0.0) {
        Nout = 0.0;
      }

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
        Pout = ((5. * P(r.ind, mesh->yend, jz) * V(r.ind, mesh->yend, jz) +
                 N(r.ind, mesh->yend, jz) * pow(V(r.ind, mesh->yend, jz), 3)) /
                    Vout -
                Nout * Vout * Vout) /
               5.;
        break;
      }
      default:
        throw BoutException("Unrecognised pressure_sheath option");
      }

      if (Pout < 0.0) {
        Pout = 0.0;
      }

      // Additional cooling
      BoutReal q;

      if (sheath_outflow) {
        // Fluid equations provide some heat flux
        q = (sheath_gamma - 6) * T(r.ind, mesh->yend, jz) * flux;
      } else {
        q = sheath_gamma * Nout * T(r.ind, mesh->yend, jz) *
            sqrt(T(r.ind, mesh->yend, jz) / AA);
      }
        
      // Multiply by cell area to get power
      BoutReal heatflux =
          q * (coord->J(r.ind, mesh->yend) + coord->J(r.ind, mesh->yend + 1)) /
          (sqrt(coord->g_22(r.ind, mesh->yend)) +
           sqrt(coord->g_22(r.ind, mesh->yend + 1)));

      // Divide by volume of cell, and 2/3 to get pressure
      ddt(P)(r.ind, mesh->yend, jz) -=
        0.5 * /// <- Because Te = Ti and this is ion only
          (2. / 3) * heatflux /
          (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));

      // Set boundary half-way between cells
      for (int jy = mesh->yend + 1; jy < mesh->LocalNy; jy++) {

        // velocity
        V(r.ind, jy, jz) = 2. * Vout - V(r.ind, mesh->yend, jz);

        // density
        N(r.ind, jy, jz) = 2 * Nout - N(r.ind, mesh->yend, jz);

        if (sheath_outflow) {
          // NV. This can be negative, so set this to the flux
          // going out of the domain (zero gradient)
          // NOTE: NV includes the ion mass AA
          NV(r.ind, jy, jz) = AA * Nout * Vout;
        } else {
          NV(r.ind, jy, jz) = 2. * AA * Nout * Vout - NV(r.ind, mesh->yend, jz);
        }

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

  /////////////////////////////////////////////////////
  // Sources
  if (density_upstream > 0.0) {

    TRACE("Density upstream");

    BoutReal source;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      int jz = 0;

      // Density source, so dn/dt = source
      BoutReal error = density_upstream - N(r.ind, mesh->ystart, jz);

      ASSERT2(finite(error));
      ASSERT2(finite(density_error_integral));

      // PI controller, using crude integral of the error
      if (density_error_lasttime < 0.0) {
        // First time
        density_error_lasttime = time;
        density_error_last = error;
      }

      // Integrate using Trapezium rule
      if (time > density_error_lasttime) { // Since time can decrease
        density_error_integral += (time - density_error_lasttime) * 0.5 *
                                  (error + density_error_last);
      }

      if ((density_error_integral < 0.0) && density_integral_positive) {
        // Limit density_error_integral to be >= 0
        density_error_integral = 0.0;
      }

      // Calculate source from combination of error and integral
      source = density_controller_p * error +
               density_controller_i * density_error_integral;

      density_error_last = error;
      density_error_lasttime = time;
    }

    if ((source < 0.0) && density_source_positive) {
      source = 0.0; // Don't remove particles
    }

    // Broadcast the value of source from processor 0
    MPI_Bcast(&source, 1, MPI_DOUBLE, 0, BoutComm::get());
    ASSERT2(finite(source));

    // Scale NeSource
    NeSource = source * NeSource0;
  }

  /////////////////////////////////////////////////////
  // Density

  // NOTE: The factor of 2 comes from the electrons
  //       Probaby need some other way to calculate maximum wave speed
  Field3D a = sqrt(gamma_sound * 2. * T / AA); // Local sound speed

  {
    TRACE("Density");

    ddt(N) = -Div_par_FV_FS(N, V, a, bndry_flux_fix);

    ddt(N) += NeSource; // External volume source

    if (anomalous_D > 0.0) {
      ddt(N) += Div_par_diffusion(anomalous_D, N);
    }
  }

  /////////////////////////////////////////////////////
  // Momentum
  
  {
    TRACE("Momentum");
    ddt(NV) =
      - Div_par_FV_FS(NV, V, a, bndry_flux_fix) // Momentum flow
      - Grad_par(P)
      ;
    
    if (viscos > 0.) {
      ddt(NV) += viscos * Div_par_diffusion_index(V);
    }

    if (anomalous_D > 0.0) {
      ddt(NV) += Div_par_diffusion(anomalous_D * V, N);
    }

    if (ion_viscosity) {
      // Braginskii ion viscosity
      BoutReal mi_me = SI::Mp / SI::Me;
      Field3D tau_i = sqrt(2 * mi_me) * Omega_ci * tau_e0 * pow(T, 1.5) / N;
      eta_i = (4. / 3) * 0.96 * N * tau_i * T;
      eta_i.applyBoundary("neumann");

      ddt(NV) += Div_par_diffusion(eta_i, V);
    }
  }

  /////////////////////////////////////////////////////
  // Pressure

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
  for (auto i : ddt(N).getRegion(RGN_NOBNDRY)) {
    if ((N[i] < 1e-5) && (ddt(N)[i] < 0.0)) {
      ddt(N)[i] = 0.0;
      ddt(NV)[i] = 0.0;
      ddt(P)[i] = 0.0;
    }
  }
}
