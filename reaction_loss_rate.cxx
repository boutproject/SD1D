// Hydrogen charge exchange

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "unused.hxx"

#include "reaction.hxx"

class ReactionLossRate : public Reaction {
public:
  ReactionLossRate(Options *options) {
    AUTO_TRACE();
    name = (*options)["name"].as<std::string>();
    
    // Loss rate in units [1/s]
    loss_rate = (*options)["loss_rate"].withDefault(0.0);
  }

  void updateSpecies(const SpeciesMap &species, BoutReal UNUSED(Tnorm), BoutReal UNUSED(Nnorm),
                     BoutReal UNUSED(Cs0), BoutReal Omega_ci) {

    AUTO_TRACE();

    // Normalised loss rate
    BoutReal loss_rate_norm = loss_rate / Omega_ci;
    
    // Get the species
    auto &atoms = *species.at(name);
    
    // Extract required variables
    Field3D N{atoms.N}, NV{atoms.NV}, P{atoms.P};
    
    Coordinates *coord = mesh->getCoordinates();

    Sn = 0.0;
    Snv = 0.0;
    Se = 0.0;
    
    CELL_AVERAGE(i,                          // Index variable
                 Sn.getRegion(RGN_NOBNDRY),  // Index and region (input)
                 coord,                      // Coordinate system (input)
                 weight,                     // Quadrature weight variable
                 N, NV, P) {   // Field variables

      // Loss of particles
      Sn[i] -= weight * loss_rate_norm * N;

      // Loss of momentum
      Snv[i] -= weight * loss_rate_norm * NV;

      // Loss of internal energy
      Se[i] -= weight * loss_rate_norm * (3./2) * P;
    }
  }
  
  SourceMap densitySources() override {
    return {{name, Sn}};
  }
  SourceMap momentumSources() {
    return {{name, Snv}};
  }
  SourceMap energySources() {
    return {{name, Se}};
  }

  std::string str() const { return "Neutral loss"; }

private:
  std::string name; // Species name e.g. "h" or "h+"
  
  Field3D Sn, Snv, Se; // Sources of particles, momentum, internal energy
  
  BoutReal loss_rate; // Neutral loss rate [1/s]
};

namespace {
RegisterInFactory<Reaction, ReactionLossRate> register_nl("fixed_loss_rate");
}
