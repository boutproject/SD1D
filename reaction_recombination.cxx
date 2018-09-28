#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "globals.hxx"
#include "unused.hxx"
#include "output.hxx"

#include "reaction.hxx"
#include "radiation.hxx"

class ReactionHydrogenRecombination : public Reaction {
public:
  ReactionHydrogenRecombination(Options *options) {
    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT4(Rrec, Erec, Frec, Srec);
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm,
                     BoutReal Nnorm, BoutReal UNUSED(Cs0), BoutReal Omega_ci) {
    TRACE("ReactionHydrogenRecombination::updateSpecies");
    
    // Get the species
    Field3D Ti, Vi;
    try {
      auto &ions = species.at("h+");
      Ti = ions.T;
      Vi = ions.V;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'h+' species");
    }
    
    // Extract required variables
    Field3D Ne, Te;
    try {
      auto &elec = species.at("e");
      Ne = elec.N;
      Te = elec.T;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'e' species");
    }
    
    Coordinates *coord = mesh->coordinates();
    
    Rrec = 0.0;
    Erec = 0.0;
    Frec = 0.0;
    Srec = 0.0;
    
    CELL_AVERAGE(i,                        // Index variable
                 Rrec.region(RGN_NOBNDRY),  // Index and region (input)
                 coord,                    // Coordinate system (input)
                 weight,                   // Quadrature weight variable
                 Ne, Te, Ti, Vi) { // Field variables

      BoutReal R = hydrogen.recombination(Ne * Nnorm, Te * Tnorm) *
        SQ(Ne) * Nnorm / Omega_ci;
      
      // Rrec is radiated energy.
      // Factor of 1.09 so that recombination becomes an energy source
      // at 5.25eV
      Rrec[i] += weight * (1.09 * Te - 13.6 / Tnorm) * R;
      
      // Erec is energy transferred to neutrals  
      Erec[i] += weight * (3. / 2) * Ti * R;

      // Momentum transferred to neutrals
      Frec[i] += weight * Vi * R;

      // Particles transferred to neutrals
      Srec[i] += weight * R;
    }
  }
  
  SourceMap densitySources() override {
    return {{"h+", -Srec},  // Sink of hydrogen ions
            {"h", Srec}};   // Source of hydrogen atoms
  }
  SourceMap momentumSources() {
    return {{"h+", -Frec}, // Plasma ions
            {"h", Frec}};  // Neutral hydrogen atoms
  }
  SourceMap energySources() {
    return {{"h+", -Erec}, // Plasma ion energy transferred to neutrals
            {"e", -Rrec},  // Electron energy radiated
            {"h", Erec}};  // Neutral hydrogen atoms
  }

  std::string str() const { return "Hydrogen recombination"; }
private:
  Field3D Rrec, Erec, Frec, Srec;

  UpdatedRadiatedPower hydrogen; // Atomic rates
};

namespace {
RegisterInFactory<Reaction, ReactionHydrogenRecombination> register_rc("recombination");
}
