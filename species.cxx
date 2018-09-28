
#include "species.hxx"
#include "div_ops.hxx"

FluidSpecies::FluidSpecies(std::string name,
                           Options *opt,
                           Solver *solver,
                           bool restarting) : name(std::move(name)) {
  
  OPTION(opt, gamma_sound, 5. / 3); // Ratio of specific heats
  OPTION(opt, bndry_flux_fix, false);

  solver->add(N, ("N"+name).c_str());
  solver->add(P, ("P"+name).c_str());
  solver->add(NV, ("NV"+name).c_str());
}

void FluidSpecies::evolve(BoutReal UNUSED(time)) {
  TRACE("Species::evolve");

  // Communicate evolving variables
  mesh->communicate(N, NV, P);

  // Floor small values
  P = floor(P, 1e-10);
  N = floor(N, 1e-10);

  Field3D Nlim = floor(N, 1e-5);

  V = NV / N; // Velocity
  T = P / N;  // Temperature

  Field3D a = sqrt(gamma_sound * T); // Local sound speed

  {
    TRACE("Density");

    ddt(N) = -Div_par_FV_FS(N, V, a, bndry_flux_fix);
  }
  {
    TRACE("Momentum");
    ddt(NV) = -Div_par_FV_FS(NV, V, a, bndry_flux_fix) // Momentum flow
               - Grad_par(P);
  }
  {
    TRACE("Pressure");

    ddt(P) += -Div_par_FV_FS(P, V, a, bndry_flux_fix) // Advection
              - (2. / 3) * P * Div_par(V)             // Compression
        ;
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
