#include "reaction.hxx"
#include "utils.hxx"
#include "unused.hxx"

/// Carbon in coronal equilibrium
/// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
class ReactionHutchinsonCarbon : public Reaction {
public:
  ReactionHutchinsonCarbon(Options* UNUSED(options)) {}
  
  void updateSpecies(const SpeciesMap &species, BoutReal Tn, BoutReal Nn,
                     BoutReal UNUSED(Vn), BoutReal UNUSED(Freq)) {
    TRACE("ReactionHutchinsonCarbon::updateSpecies");
    
    // Electron temperature and density
    try {
      Te = species.at("e").T;
      ne = species.at("e").N;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'e' species");
    }
    
    // Carbon density
    try {
      ni = species.at("c").N;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'c' species");
    }

    // Store normalisations
    Tnorm = Tn;
    Nnorm = Nn;
    // Velocity and frequency normalisation not needed
  }

  SourceMap energySources() {
    // Calculate radiated power
    // Note Te, ne, ni normalised
    Field3D radiation = SQ(Nnorm) * ne * ni * 2e-31 * pow(Te * Tnorm / 10., 3) /
                        (1. + pow(Te * Tnorm / 10., 4.5));

    // Loss of energy for electrons
    return {{"e", -radiation}};
  }

  std::string str() const { return "Coronal carbon radiation (Hutchinson model)"; }
  
private:
  BoutReal Nnorm, Tnorm; ///< Density, temperature normalisations
  Field3D Te, ne, ni;
};

namespace {
  RegisterInFactory<Reaction, ReactionHutchinsonCarbon> register_ch("c_hutchinson");
}

