#include "fluid_full.hxx"

#include "div_ops.hxx"

void FluidFull::init(Solver *solver,
                     bool restarting,
                     const string &suffix, 
                     BoutReal Tnorm, BoutReal Nnorm,
                     BoutReal Cs0, BoutReal rho_s0) : {

  // Specify evolving variables. Each fluid must have a distinct suffix to avoid
  // name clashes
  solver->add(N, string("N")+suffix);
  solver->add(P, string("P")+suffix);
  solver->add(NVi, string("NV")+suffix);
  
  
}

/// Calculate time derivatives
void FluidFull::rhs(BoutReal time) {
  mesh->communicate(N, NV, P);
  
  // Floor small values
  P = floor(P, 1e-12);
  N = floor(N, 1e-10);

  V = NV / N;
  T = P / N;
  
  // Apply boundary conditions
  
  
  // Calculate local sound speed
  // This is used in upwinding schemes to provide some dissipation
  Field3D cs = sqrt( gamma*P/N ); // Local sound speed
  
  // Density equation
  ddt(N) = - Div_par_FV_FS(N, V, cs, true);

  // Parallel momentum
  ddt(NV) =
    - Div_par_FV_FS(NVi, V, cs, false) // Momentum flow
    - Grad_par(P) / AA; // Pressure force
  ;
  
  // Pressure
  
  ddt(P) += 
    - Div_par_FV_FS(P, V, cs, true)       // Advection
    - P*Div_par(V)*(gamma-1.0)            // Compression
    ;

  if (diffusion_multiplier > 0.0) {
    // Include a parallel diffusion driven by pressure
    
  }
}

BoutReal FluidFull::density(const DataIterator &i) {
  return N[i];
}

BoutReal FluidFull::temperature(const DataIterator &i) {
  return P[i] / N[i];
}

BoutReal FluidFull::velocity(const DataIterator &i) {
  return NV[i] / N[i];
}

void FluidFull::addParticles(const DataIterator &i, BoutReal Sn) {
  ddt(N)[i] += Sn;
}

void FluidFull::addEnergy(const DataIterator &i, BoutReal Se) {
  ddt(P)[i] += Se/gamma;
}

void FluidFull::addMomentum(const DataIterator &i, BoutReal Sm) {
  ddt(NV) += Sm/AA;
}
