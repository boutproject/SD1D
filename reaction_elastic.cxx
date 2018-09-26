/////////////////////////////////////////////////////////
// Ion-neutral elastic scattering
//
// Post "A Review of Recent Developments in Atomic Processes for
// Divertors and Edge Plasmas" PSI review paper
//       https://arxiv.org/pdf/plasm-ph/9506003.pdf
// Relative velocity of two particles in a gas
// is sqrt(8kT/pi mu) where mu = m_A*m_B/(m_A+m_B)
// here ions and neutrals have same mass,
// and the ion temperature is used

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "globals.hxx"
#include "reaction.hxx"
#include "unused.hxx"

#include "output.hxx"

class ReactionElasticScattering : public Reaction {
public:
  ReactionElasticScattering() { SAVE_REPEAT2(Fel, Eel); }

  void updateSpecies(const SpeciesMap &species, BoutReal UNUSED(Tnorm),
                     BoutReal Nnorm, BoutReal Cs0, BoutReal Omega_ci) {

    TRACE("ReactionElasticScattering::updateSpecies");

    // Plasma ions
    Field3D Ti, Ni, Vi;
    try {
      const auto &ions = species.at("h+");
      
      Ti = ions.T;
      Ni = ions.N;
      Vi = ions.V;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'h+' species");
    }

    // Neutral atoms
    Field3D Tn, Nn, Vn;
    try {
      const auto &atoms = species.at("h");
      Tn = atoms.T;
      Nn = atoms.N;
      Vn = atoms.V;
    } catch (const std::out_of_range &e) {
      throw BoutException("No 'h' species");
    }
    
    Coordinates *coord = mesh->coordinates();

    Fel = 0.0;
    Eel = 0.0;
    
    CELL_AVERAGE(i,                        // Index variable
                 Fel.region(RGN_NOBNDRY),  // Index and region (input)
                 coord,                    // Coordinate system (input)
                 weight,                   // Quadrature weight variable
                 Ti, Ni, Vi, Tn, Nn, Vn) { // Field variables
      
      // Rate (normalised)
      BoutReal R =
          a0 * Ni * Nn * Cs0 * sqrt((16. / PI) * Ti) * Nnorm / Omega_ci;

      // Elastic transfer of momentum
      Fel[i] += weight * (Vi - Vn) * R;

      // Elastic transfer of thermal energy
      Eel[i] += weight * (3. / 2) * (Ti - Tn) * R;
    }
  }

  SourceMap momentumSources() {
    return {{"h+", -Fel}, // Deuterium (plasma ions)
            {"h", Fel}}; // Neutral atoms
  }
  SourceMap energySources() {
    return {{"h+", -Eel}, // Deuterium (plasma ions)
            {"h", Eel}}; // Neutral atoms
  }

  std::string str() const { return "Ion-neutral elastic scattering"; }

private:
  BoutReal a0 = 3e-19; // Effective cross-section [m^2]

  Field3D Fel, Eel;
};

namespace {
RegisterInFactory<Reaction, ReactionElasticScattering> register_el("elastic");
}
