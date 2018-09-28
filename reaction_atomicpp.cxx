
// OpenADAS interface Atomicpp by T.Body
#include "atomicpp/ImpuritySpecies.hxx"
#include "atomicpp/Prad.hxx"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "unused.hxx"

#include "reaction.hxx"

// Use OpenADAS data through Atomicpp
class ReactionAtomicppCoronal : public Reaction {
public:
  ReactionAtomicppCoronal(Options *options) {

    // Find out which species to model
    OPTION(options, impurity_species, "c");
    impurity = new ImpuritySpecies(impurity_species);

    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT(Rzrad);
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm, BoutReal Nnorm,
                     BoutReal UNUSED(Vn), BoutReal UNUSED(Omega_ci)) {
    TRACE("ReactionHutchinsonCarbon::updateSpecies");

    Field3D Te, Ne, Nz, Nn;
    // Electron temperature and density
    try {
      Te = species.at("e").T;
      Ne = species.at("e").N;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'e' species");
    }

    // Neutral atom density
    try {
      Nn = species.at("h").N;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'h' species");
    }

    // Impurity density
    try {
      Nz = species.at(impurity_species).N;
    } catch (const std::out_of_range &e) {
      throw BoutException("No '%s' species", impurity_species.c_str());
    }

    Coordinates *coord = mesh->coordinates();

    Rzrad = 0.0;

    CELL_AVERAGE(i,                         // Index variable
                 Rzrad.region(RGN_NOBNDRY), // Index and region (input)
                 coord,                     // Coordinate system (input)
                 weight,                    // Quadrature weight variable
                 Ne, Nz, Nn, Te) {          // Field variables
      Rzrad[i] +=
          weight * computeRadiatedPower(*impurity,
                                        Te * Tnorm, // electron temperature [eV]
                                        Ne * Nnorm, // electron density [m^-3]
                                        Nz * Nnorm, // impurity density [m^-3]
                                        Nn * Nnorm); // Neutral density [m^-3]
    }
  }

  SourceMap energySources() {
    return {{"e", -Rzrad}}; // Electron energy sink
  }

  std::string str() const {
    return "Coronal radiation from species '" + impurity_species +
           "' (Atomic++ OpenADAS data)";
  }

private:
  Field3D Rzrad;
  string impurity_species;

  ImpuritySpecies *impurity; // Atomicpp impurity
};

namespace {
RegisterInFactory<Reaction, ReactionAtomicppCoronal>
    register_at("atomic++coronal");
}
