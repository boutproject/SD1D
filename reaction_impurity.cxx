#include "bout/constants.hxx"
#include "utils.hxx"
#include "unused.hxx"

#include "reaction.hxx"

/// Carbon in coronal equilibrium
/// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
class ReactionHutchinsonCarbon : public Reaction {
public:
  ReactionHutchinsonCarbon(Options* UNUSED(options)) {}
  
  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm, BoutReal Nnorm,
                     BoutReal UNUSED(Vn), BoutReal Omega_ci) {
    TRACE("ReactionHutchinsonCarbon::updateSpecies");
    
    // Electron temperature and density
    Field3D Te, ne;
    try {
      Te = species.at("e")->T;
      ne = species.at("e")->N;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'e' species");
    }
    
    // Carbon density
    Field3D ni;
    try {
      ni = species.at("c")->N;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'c' species");
    }

    // Calculate radiated power
    // Note Te, ne, ni normalised
    radiation = SQ(Nnorm) * ne * ni * 2e-31 * pow(Te * Tnorm / 10., 3) /
      (1. + pow(Te * Tnorm / 10., 4.5))
      / (SI::qe * Tnorm * Nnorm * Omega_ci); // Normalise
    
  }

  SourceMap energySources() {
    // Loss of energy for electrons
    return {{"e", -radiation}};
  }

  std::string str() const { return "Coronal carbon radiation (Hutchinson model)"; }
  
private:
  Field3D radiation;
};

namespace {
  RegisterInFactory<Reaction, ReactionHutchinsonCarbon> register_ch("c_hutchinson");
}

