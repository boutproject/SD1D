// Hydrogen charge exchange

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/fv_ops.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "unused.hxx"

#include "radiation.hxx"
#include "reaction.hxx"

using bout::globals::mesh;

class ReactionNeutralDiffusion : public Reaction {
public:
  ReactionNeutralDiffusion(Options *options) {
    AUTO_TRACE();
    
    bool diagnose = (*options)["diagnose"].withDefault(false);
    if (diagnose) {
      SAVE_REPEAT(Dn);
    }

    // Floor on neutral temperature for diffusion calculation
    // Set to approximate Frank-Condon energy [eV]
    tn_floor = (*options)["tn_floor"].withDefault(3.5);
    
    // Maximum mean free path [m]
    max_mfp = (*options)["max_mfp"].withDefault(0.1);

    // Scale neutral gas diffusion
    dneut = (*options)["dneut"].withDefault(1.0); 
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm, BoutReal Nnorm,
                     BoutReal Cs0, BoutReal Omega_ci) {

    BoutReal rho_s0 = Cs0 / Omega_ci;
    
    // Get the species
    auto &atoms = *species.at("h");
    auto &ions = *species.at("h+");
    auto &electrons = *species.at("e");
    
    // Extract required variables
    Field3D Ni{ions.N}, Ti{ions.T};
    Field3D Nn{atoms.N}, Tn{atoms.T}, NVn{atoms.NV}, Vn{atoms.V}, Pn{atoms.P};
    Field3D Ne{electrons.N}, Te{electrons.T};
    
    Coordinates *coord = mesh->getCoordinates();

    Dn = 0.0;
    
    CELL_AVERAGE(i,                          // Index variable
                 Dn.getRegion(RGN_NOBNDRY),  // Index and region (input)
                 coord,                      // Coordinate system (input)
                 weight,                     // Quadrature weight variable
                 Ni, Nn, Ne, Ti, Te, Tn) {   // Field variables

      // Charge exchange frequency, normalised to ion cyclotron
      // frequency
      BoutReal sigma_cx =
        Ni * Nnorm * hydrogen.chargeExchange(Ti * Tnorm) / Omega_ci;

      // Ionisation frequency
      BoutReal sigma_iz =
        Ne * Nnorm * hydrogen.ionisation(Te * Tnorm) / Omega_ci;

      // Neutral thermal velocity
      BoutReal tn = Tn;
      if (tn < tn_floor / Tnorm) {
        tn = tn_floor / Tnorm;
      }
      BoutReal vth_n = sqrt(tn); // Normalised to Cs0
      
      // Neutral-neutral mean free path
      BoutReal a0 = PI * SQ(5.29e-11);
      BoutReal lambda_nn = 1. / (Nnorm * Nn * a0); // meters
      if (lambda_nn > max_mfp) {
        // Limit maximum mean free path
        lambda_nn = max_mfp;
      }
      
      lambda_nn /= rho_s0; // Normalised length to rho_s0
      // Neutral-Neutral collision rate
      BoutReal sigma_nn = vth_n / lambda_nn;
      
      // Total neutral collision frequency
      BoutReal sigma = sigma_cx + sigma_iz + sigma_nn;
      
      // Neutral gas diffusion coefficient
      if (dneut > 0.0) {
        BoutReal update = weight * dneut * SQ(vth_n) / sigma;
        
        Dn[i] += update;
        kappa_n[i] += update * Nn; // heat conduction
      }
    }

    // Apply boundaries, communicate
    kappa_n.applyBoundary("neumann");
    Dn.applyBoundary("dirichlet_o2");
    mesh->communicate(kappa_n, Dn);

    // Relationship between heat conduction and viscosity for neutral
    // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
    // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
    // Transport Processes in Gases", 1972
    //
    eta_n = (2. / 5) * kappa_n;

    // Cross-field diffusion calculated from pressure gradient
    Field3D logPn = log(floor(Pn, 1e-7));
    logPn.applyBoundary("neumann");
    
    // Particle diffusion
    Sn = FV::Div_par_K_Grad_par(Dn * Nn, logPn);

    // Momentum diffusion
    Snv = FV::Div_par_K_Grad_par(NVn * Dn, logPn)
      + FV::Div_par_K_Grad_par(eta_n, Vn);

    // Heat conduction
    Se = FV::Div_par_K_Grad_par(kappa_n, Tn) // Parallel
      + FV::Div_par_K_Grad_par(Dn * (3./2)*Pn, logPn); // Perpendicular diffusion
    
  }
  
  SourceMap densitySources() override {
    return {{"h", Sn}};
  }
  SourceMap momentumSources() {
    return {{"h", Snv}};
  }
  SourceMap energySources() {
    return {{"h", Se}};
  }

  std::string str() const { return "Neutral diffusion"; }

private:
  Field3D Dn, kappa_n, eta_n; // Particle diffusion, thermal conduction, viscosity
  
  Field3D Sn, Snv, Se; // Sources of particles, momentum, internal energy
  
  BoutReal tn_floor; // Floor on neutral temperature [eV]
  BoutReal max_mfp; // Maximum mean free path [m]
  BoutReal dneut; // Multiply diffusion rate
  
  UpdatedRadiatedPower hydrogen; // Atomic rates
};

namespace {
RegisterInFactory<Reaction, ReactionNeutralDiffusion> register_nd("neutral_diffusion");
}
