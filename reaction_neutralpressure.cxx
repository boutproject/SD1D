#include "reaction.hxx"
#include "difops.hxx"
#include "options.hxx"

/// A simple model for friction with neutrals
/// Sets the force on the plasma to be the gradient of the neutral pressure
class ReactionNeutralPressureForce : public Reaction {
public:
  ReactionNeutralPressureForce(Options *UNUSED(options)) {}

  void updateSpecies(const SpeciesMap &species, BoutReal UNUSED(Tnorm),
                     BoutReal UNUSED(Nnorm), BoutReal UNUSED(Cs0), BoutReal UNUSED(Omega_ci)) {
    TRACE("ReactionNeutralPressureForce::updateSpecies");
    
    Field3D Pn;
    try {
      Pn = species.at("h").P; // Neutral pressure
    } catch(const std::out_of_range &e) {
      throw BoutException("No 'h' species");
    }
    Fpn = Grad_par(Pn);
  }
  SourceMap momentumSources() {
    return {{"h+", -Fpn}, // Plasma ions
            {"h", Fpn}};  // Neutral hydrogen atoms
  }
private:
  Field3D Fpn;
};

namespace {
RegisterInFactory<Reaction, ReactionNeutralPressureForce> register_rc("neutralpressureforce");
}
