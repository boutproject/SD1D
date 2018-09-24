#include "reaction.hxx"
#include "utils.hxx"

/// Carbon in coronal equilibrium
/// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
class ReactionHutchinsonCarbon : public Reaction {
public:
  void updateSpecies(const SpeciesMap &species, BoutReal Tn, BoutReal Nn,
                     BoutReal Vn) {
    // Electron temperature and density
    Te = species.at("e").T;
    ne = species.at("e").N;
    // Carbon density
    ni = species.at("c").N;

    // Store normalisations
    Tnorm = Tn;
    Nnorm = Nn;
    // Velocity normalisation not needed
  }

  SourceMap energySources() {
    // Calculate radiated power
    // Note Te, ne, ni normalised
    Field3D radiation = SQ(Nnorm) * ne * ni * 2e-31 * pow(Te * Tnorm / 10., 3) /
                        (1. + pow(Te * Tnorm / 10., 4.5));

    // Loss of energy for electrons
    return {{"e", -radiation}};
  }

private:
  BoutReal Nnorm, Tnorm; ///< Density, temperature normalisations
  Field3D Te, ne, ni;
};

namespace {
  RegisterInFactory<Reaction, ReactionHutchinsonCarbon> register_ch("c_hutchinson");
}

