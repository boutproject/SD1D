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
    Field3D Ti = species.at("i").T;
    Field3D Ni = species.at("i").N;
    Field3D Vi = species.at("i").V;

    // Neutral atoms
    Field3D Tn = species.at("n").T;
    Field3D Nn = species.at("n").N;
    Field3D Vn = species.at("n").V;

    Coordinates *coord = mesh->coordinates();

    Fel = 0.0;
    Eel = 0.0;
    
    INTEGRATE(i,                        // Index variable
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
    return {{"i", -Fel}, // Deuterium (plasma ions)
            {"n", Fel}}; // Neutral atoms
  }
  SourceMap energySources() {
    return {{"i", -Eel}, // Deuterium (plasma ions)
            {"n", Eel}}; // Neutral atoms
  }

  std::string str() const { return "Ion-neutral elastic scattering"; }

private:
  BoutReal a0 = 3e-19; // Effective cross-section [m^2]

  Field3D Fel, Eel;
};

namespace {
RegisterInFactory<Reaction, ReactionElasticScattering> register_el("elastic");
}
