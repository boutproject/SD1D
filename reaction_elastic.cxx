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

    Fel.allocate();
    Eel.allocate();

    for (auto i : Fel.region(RGN_NOBNDRY)) {

      InterpCell cell_Ti(Ti, i);
      InterpCell cell_Ni(Ni, i);
      InterpCell cell_Vi(Vi, i);

      InterpCell cell_Tn(Tn, i);
      InterpCell cell_Nn(Nn, i);
      InterpCell cell_Vn(Vn, i);

      // Jacobian (Cross-sectional area)
      InterpCell cell_J(coord->J, i);

      // Rates (normalised)
      BoutReal R_el_L = a0 * cell_Ni.l * cell_Nn.l * Cs0 *
                        sqrt((16. / PI) * cell_Ti.l) * Nnorm / Omega_ci;
      BoutReal R_el_C = a0 * cell_Ni.c * cell_Nn.c * Cs0 *
                        sqrt((16. / PI) * cell_Ti.c) * Nnorm / Omega_ci;
      BoutReal R_el_R = a0 * cell_Ni.r * cell_Nn.r * Cs0 *
                        sqrt((16. / PI) * cell_Ti.r) * Nnorm / Omega_ci;

      // Elastic transfer of momentum
      Fel[i] = (cell_J.l * (cell_Vi.l - cell_Vn.l) * R_el_L +
                4. * cell_J.c * (cell_Vi.c - cell_Vn.c) * R_el_C +
                cell_J.r * (cell_Vi.r - cell_Vn.r) * R_el_R) /
               (6. * cell_J.c);

      // Elastic transfer of thermal energy
      Eel[i] = (3. / 2) *
               (cell_J.l * (cell_Ti.l - cell_Tn.l) * R_el_L +
                4. * cell_J.c * (cell_Ti.l - cell_Tn.l) * R_el_C +
                cell_J.r * (cell_Ti.r - cell_Tn.r) * R_el_R) /
               (6. * cell_J.c);
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
