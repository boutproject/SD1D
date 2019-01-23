// Hydrogen charge exchange

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "unused.hxx"

#include "radiation.hxx"
#include "reaction.hxx"

class ReactionHydrogenCX : public Reaction {
public:
  ReactionHydrogenCX(Options *options) {
    TRACE("ReactionHydrogenCX(Options*)");

    OPTION(options, charge_exchange_escape, false);
    if (charge_exchange_escape) {
      // Charge exchanged neutrals escape the plasma and are
      // redistributed.

      // Fraction of energy which returns
      OPTION(options, charge_exchange_return_fE, 1.0);

      // Calculate neutral gas redistribution weights over the domain
      std::string redist_string;
      options->get("redist_weight", redist_string, "1.0");

      FieldFactory ffact(mesh);
      redist_weight = ffact.create2D(redist_string, options);
      BoutReal localweight = 0.0;
      Coordinates *coord = mesh->getCoordinates();
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        localweight += redist_weight(mesh->xstart, j) *
                       coord->J(mesh->xstart, j) * coord->dy(mesh->xstart, j);
      }

      ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator

      // Calculate total weight by summing over all processors
      BoutReal totalweight;
      MPI_Allreduce(&localweight, &totalweight, 1, MPI_DOUBLE, MPI_SUM, ycomm);
      // Normalise redist_weight so sum over domain:
      //
      // sum ( redist_weight * J * dy ) = 1
      //
      redist_weight /= totalweight;
    }

    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT2(Fcx, Ecx);
      if (charge_exchange_escape) {
        SAVE_REPEAT(Scx);
      }
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm, BoutReal Nnorm,
                     BoutReal UNUSED(Cs0), BoutReal Omega_ci) {

    // Get the species
    auto &atoms = *species.at("h");
    auto &ions = *species.at("h+");

    // Extract required variables
    Field3D Ni{ions.N}, Ti{ions.T}, Vi{ions.V};
    Field3D Nn{atoms.N}, Tn{atoms.T}, Vn{atoms.V};

    Coordinates *coord = mesh->getCoordinates();

    Ecx = 0.0;
    Fcx = 0.0;
    if (charge_exchange_escape) {
      Scx = 0.0;
      Scx_E = 0.0;
    }

    CELL_AVERAGE(i,                          // Index variable
                 Ecx.getRegion(RGN_NOBNDRY), // Index and region (input)
                 coord,                      // Coordinate system (input)
                 weight,                     // Quadrature weight variable
                 Ni, Nn, Ti, Tn, Vi, Vn) {   // Field variables

      BoutReal R =
          Ni * Nn * hydrogen.chargeExchange(Ti * Tnorm) * (Nnorm / Omega_ci);

      // Ecx is energy transferred from ions to neutrals
      Ecx[i] += weight * (3. / 2) * (Ti - Tn) * R;

      // Fcx is friction between plasma and neutrals
      Fcx[i] += weight * (Vi - Vn) * R;

      if (charge_exchange_escape) {
        // Scx is a redistribution of fast neutrals due to charge exchange
        // Acts as a sink of plasma density
        Scx[i] -= weight * R;

        // Energy lost from the neutrals
        // Note: Has ion temperature since it's a charge-exchanged neutral
        Scx_E[i] -= weight * (3. / 2) * Ti * R;
      }
    }

    if (charge_exchange_escape) {
      // Fast CX neutrals lost from plasma.
      // These are redistributed, along with a fraction of their energy

      BoutReal Scx_Ntot = 0.0;
      BoutReal Scx_Etot = 0.0;
      for (auto &i : Scx.getRegion(RGN_NOBNDRY)) {
        Scx_Ntot += Scx[i] * coord->J[i] * coord->dy[i];
        Scx_Etot += Scx_E[i] * coord->J[i] * coord->dy[i];
      }

      // Now sum on all processors
      BoutReal send[2] = {Scx_Ntot, Scx_Etot};
      BoutReal recv[2];
      MPI_Allreduce(send, recv, 2, MPI_DOUBLE, MPI_SUM, ycomm);
      Scx_Ntot = recv[0];
      Scx_Etot = recv[1];

      // Scale the energy of the returning CX neutrals
      Scx_Etot *= charge_exchange_return_fE;

      // Use the normalised redistribuion weight
      // sum ( redist_weight * J * dy ) = 1
      for (auto &i : Scx.getRegion(RGN_NOBNDRY)) {
        // Note: Scx_Ntot, Scx_Etot < 0 here
        Scx[i] -= Scx_Ntot * redist_weight[i];
        Scx_E[i] -= Scx_Etot * redist_weight[i];
      }
    }
  }
  SourceMap densitySources() override {
    if (charge_exchange_escape) {
      return {{"h", Scx}}; // Redistribute neutrals
    }
    return {};
  }
  SourceMap momentumSources() {
    if (charge_exchange_escape) {
      return {{"h+", -Fcx}}; // Plasma ions
      // Note: neutrals escape, losing momentum to the wall
    }
    return {{"h+", -Fcx}, // Plasma ions
            {"h", Fcx}};  // Neutral hydrogen atoms
  }
  SourceMap energySources() {
    if (charge_exchange_escape) {
      return {{"h+", -Ecx},        // Plasma ion energy transferred to neutrals
              {"h", Ecx + Scx_E}}; // Neutral hydrogen atoms
    }
    return {{"h+", -Ecx}, // Plasma ion energy transferred to neutrals
            {"h", Ecx}};  // Neutral hydrogen atoms
  }

  std::string str() const { return "Hydrogen charge exchange"; }

private:
  bool charge_exchange_escape;
  Field2D redist_weight;              // Weighting used to decide redistribution
  BoutReal charge_exchange_return_fE; // Fraction of energy carried by returning
                                      // CX neutrals

  Field3D Ecx, Fcx, Scx, Scx_E;

  UpdatedRadiatedPower hydrogen; // Atomic rates

  MPI_Comm ycomm;
};

namespace {
RegisterInFactory<Reaction, ReactionHydrogenCX> register_rc("hydrogen_cx");
}
