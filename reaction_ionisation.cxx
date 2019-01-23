#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "globals.hxx"
#include "unused.hxx"
#include "output.hxx"

#include "reaction.hxx"
#include "radiation.hxx"

class ReactionHydrogenIonisation : public Reaction {
public:
  ReactionHydrogenIonisation(Options *options) {
    TRACE("ReactionHydrogenIonisation(Options*)");

    // Ionisation energy cost
    OPTION(options, Eionize, 30);

    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT4(Riz, Eiz, Fiz, Siz);
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm,
                     BoutReal Nnorm, BoutReal UNUSED(Cs0), BoutReal Omega_ci) {
    // Get the species
    // Extract required variables
    Field3D Nn, Tn, Vn;
    try {
      auto &atoms = *species.at("h");
      Nn = atoms.N;
      Tn = atoms.T;
      Vn = atoms.V;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'h' species");
    }
    
    Field3D Ne, Te;
    try {
      auto &elec  = *species.at("e");
      Ne = elec.N;
      Te = elec.T;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'e' species");
    }
    
    Coordinates *coord = mesh->getCoordinates();

    Riz = 0.0;
    Eiz = 0.0;
    Fiz = 0.0;
    Siz = 0.0;
    
    CELL_AVERAGE(i,                        // Index variable
                 Riz.getRegion(RGN_NOBNDRY),  // Index and region (input)
                 coord,                    // Coordinate system (input)
                 weight,                   // Quadrature weight variable
                 Ne, Nn, Te, Tn, Vn) {     // Field variables

      BoutReal R =
        Ne * Nn * hydrogen.ionisation(Te * Tnorm) * Nnorm / Omega_ci;

      // Electron energy loss per ionisation
      Riz[i] += weight * (Eionize / Tnorm) * R;

      // Energy from neutral atom temperature
      Eiz[i] -= weight * (3. / 2) * Tn * R;

      // Friction due to ionisation
      Fiz[i] -= weight * Vn * R;

      // Plasma sink due to ionisation (negative)
      Siz[i] -= weight * R;
    }
  }
  
  SourceMap densitySources() override {
    return {{"h+", -Siz},  // Siz < 0 => ion source 
            {"h", Siz}};   // Siz < 0 => atom sink
  }
  SourceMap momentumSources() {
    return {{"h+", -Fiz}, // Plasma ion momentum source
            {"h", Fiz}};  // Neutral atom momentum sink
  }
  SourceMap energySources() {
    return {{"h+", -Eiz}, // Eiz < 0 => Ion energy source
            {"e", -Riz},  // Electron energy into ionisation
            {"h", Eiz}};  // Neutral atom energy sink
  }
  
  std::string str() const { return "Hydrogen ionisation"; }
  
private:
  UpdatedRadiatedPower hydrogen; // Atomic rates
  Field3D Riz, Eiz, Fiz, Siz;

  BoutReal Eionize;  // Ionisation energy loss [eV]
};

namespace {
RegisterInFactory<Reaction, ReactionHydrogenIonisation> register_iz("ionisation");
}
