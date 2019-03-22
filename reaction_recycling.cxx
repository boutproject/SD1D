// Recycling

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "unused.hxx"

#include "reaction.hxx"

using bout::globals::mesh;

class ReactionRecycling : public Reaction {
public:
  ReactionRecycling(Options *options) {
    AUTO_TRACE();
    name = (*options)["name"].withDefault<std::string>("h");

    vwall = (*options)["vwall"].withDefault(1./3); // velocity corresponding to 1/3rd Franck-Condon at wall
    frecycle = (*options)["frecycle"].withDefault(0.0); // Recycling fraction
    gaspuff = (*options)["gaspuff"].withDefault(0.0); // Additional gas flux at target
    fredistribute = (*options)["fredistribute"].withDefault(0.0); // Fraction of neutrals redistributed evenly along leg

    // Calculate the weighting for redistribution of neutrals
    
    std::string redist_string = (*options)["redist_weight"].withDefault<std::string>("1.0");

    FieldFactory ffact(mesh);
    redist_weight = ffact.create2D(redist_string, options);

    // Normalise so that the integral over the domain is 1.
    
    BoutReal localweight = 0.0;
    Coordinates *coord = mesh->getCoordinates();
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      localweight += redist_weight(mesh->xstart, j) *
        coord->J(mesh->xstart, j) * coord->dy(mesh->xstart, j);
    }

    MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator
    
    // Calculate total weight by summing over all processors
    BoutReal totalweight;
    MPI_Allreduce(&localweight, &totalweight, 1, MPI_DOUBLE, MPI_SUM, ycomm);
    // Normalise redist_weight so sum over domain:
    //
    // sum ( redist_weight * J * dy ) = 1
    //
    redist_weight /= totalweight;

    if ((*options)["diagnose"].withDefault(false)) {
      // Save additional outputs
      SAVE_REPEAT(flux_ion);
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm, BoutReal UNUSED(Nnorm),
                     BoutReal UNUSED(Cs0), BoutReal UNUSED(Omega_ci)) {

    AUTO_TRACE();
    
    // Get the species
    auto &atoms = *species.at(name);
    auto &ions = *species.at(name + "+");
    
    // Extract required variables
    Field3D Nn{atoms.N}, Vn{atoms.V};
    Field3D Ni{ions.N}, Vi{ions.V};
    
    Coordinates *coord = mesh->getCoordinates();
    
    // Sources
    SNVn = 0.0;
    SEn = 0.0;
    SNn = 0.0;
    
    BoutReal nredist;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      int jz = 0; // Z index
      int jy = mesh->yend;

      BoutReal ion_density = 0.5 * (Ni(r.ind, jy, jz) + Ni(r.ind, jy + 1, jz));
      BoutReal ion_velocity = 0.5 * (Vi(r.ind, jy, jz) + Vi(r.ind, jy + 1, jz));
      
      // Ion flux to the target
      flux_ion = ion_density * ion_velocity *
        (coord->J(r.ind, jy) + coord->J(r.ind, jy + 1)) /
        (sqrt(coord->g_22(r.ind, jy)) + sqrt(coord->g_22(r.ind, jy + 1)));

      BoutReal neutral_density = 0.5 * (Nn(r.ind, jy, jz) + Nn(r.ind, jy + 1, jz));
      BoutReal neutral_velocity = 0.5 * (Vn(r.ind, jy, jz) + Vn(r.ind, jy + 1, jz));
      
      BoutReal flux_neut = neutral_density * neutral_velocity *
        (coord->J(r.ind, jy) + coord->J(r.ind, jy + 1)) /
        (sqrt(coord->g_22(r.ind, jy)) + sqrt(coord->g_22(r.ind, jy + 1)));
      
      // Total amount of neutral gas to be added
      // This is split between flux recycled at the target, and flux
      // redistributed along the domain
      BoutReal nadd = flux_ion * frecycle + flux_neut + gaspuff;
      
      // Neutral gas arriving at the target
      BoutReal ntarget =
        (1 - fredistribute) * nadd /
        (coord->J(r.ind, mesh->yend) * coord->dy(r.ind, mesh->yend));
      
      SNn(r.ind, mesh->yend, jz) += ntarget;

      // Set velocity of neutrals coming from the wall to a fraction of
      // the Franck-Condon energy
      BoutReal Vneut = -vwall * sqrt(3.5 / Tnorm);
      SNVn(r.ind, mesh->yend, jz) += ntarget * Vneut;
      
      // Set temperature of the incoming neutrals to F-C
      SEn(r.ind, mesh->yend, jz) += (3./2) * ntarget * (3.5 / Tnorm);
      
      // Re-distribute neutrals
      nredist = fredistribute * nadd;

      // Divide flux_ion by J so that the result in the output file has
      // units of flux per m^2
      flux_ion /= coord->J(mesh->xstart, mesh->yend + 1);
    }

    if (fredistribute > 0.0) {
      // Now broadcast redistributed neutrals to other processors
      MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator
      int np;
      MPI_Comm_size(ycomm, &np); // Number of processors

      // Broadcast from final processor (presumably with target)
      // to all other processors
      MPI_Bcast(&nredist, 1, MPI_DOUBLE, np - 1, ycomm);

      // Distribute along length
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        // Neutrals into this cell
        // Note: from earlier normalisation the sum ( redist_weight * J * dy )
        // = 1 This ensures that if redist_weight is constant then the source
        // of particles per volume is also constant.
        BoutReal ncell = nredist * redist_weight(mesh->xstart, j);
        
        SNn(mesh->xstart, j, 0) += ncell;

        // No momentum

        // Set temperature of the incoming neutrals to F-C
        SEn(mesh->xstart, j, 0) += (3./2) * ncell * (3.5 / Tnorm);
      }
    }
  }
  
  SourceMap densitySources() override {
    return {{name, SNn}};
  }
  SourceMap momentumSources() {
    return {{name, SNVn}};
  }
  SourceMap energySources() {
    return {{name, SEn}};
  }

  std::string str() const { return "Neutral recycling"; }

private:
  std::string name; // Species name e.g. "h" or "h+"

  BoutReal flux_ion; // Flux of ions to target (output)

  
  BoutReal frecycle; // Recycling fraction
  BoutReal gaspuff;  // Additional source of neutral gas at the target plate
  BoutReal vwall;    // Velocity of neutrals coming from the wall
                     // as fraction of Franck-Condon energy
  BoutReal fredistribute; // Fraction of recycled neutrals re-distributed along
                          // length

  Field2D redist_weight; // Weighting used to decide redistribution
  
  Field3D SNn, SNVn, SEn; // Sources of particles, momentum, internal energy
};

namespace {
RegisterInFactory<Reaction, ReactionRecycling> register_rec("recycling");
}
