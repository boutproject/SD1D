#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "globals.hxx"
#include "unused.hxx"
#include "output.hxx"

#include "reaction.hxx"
#include "radiation.hxx"

class ReactionHydrogenExcitation : public Reaction {
public:
  ReactionHydrogenExcitation(Options *options) {
    TRACE("ReactionHydrogenExcitation(Options*)");
    
    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT(Rex);
    }
  }
  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm,
                     BoutReal Nnorm, BoutReal UNUSED(Cs0), BoutReal Omega_ci) {

    // Get the species
    auto &atoms = species.at("h");
    auto &elec  = species.at("e");

    // Extract required variables
    Field3D Ne{elec.N}, Te{elec.T};
    Field3D Nn{atoms.N};

    Coordinates *coord = mesh->coordinates();
    
    Rex = 0.0;
    
    CELL_AVERAGE(i,                        // Index variable
                 Rex.region(RGN_NOBNDRY),  // Index and region (input)
                 coord,                    // Coordinate system (input)
                 weight,                   // Quadrature weight variable
                 Ne, Nn, Te) {             // Field variables
      BoutReal R =
        Ne * Nn * hydrogen.excitation(Te * Tnorm) * Nnorm / Omega_ci / Tnorm;
      
      Rex[i] += weight * R;
    }
  }
  
  SourceMap energySources() {
    return {{"e", -Rex}}; // Electron energy sink
  }
  
  std::string str() const { return "Hydrogen excitation radiation"; }
  
private:
  UpdatedRadiatedPower hydrogen; // Atomic rates
  Field3D Rex;
  
};

namespace {
RegisterInFactory<Reaction, ReactionHydrogenExcitation> register_ex("excitation");
}
