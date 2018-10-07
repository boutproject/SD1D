///
/// The following labels are used for common species:
///
///   e    Electrons
///   i    Ions (main species; Deuterium)
///   n    Neutral atoms (main species; Deuterium)
///   c    Carbon
///

#ifndef __SPECIES_HXX__
#define __SPECIES_HXX__

#include "bout/generic_factory.hxx"

#include "bout/solver.hxx"

#include <map>
#include <string>

/// Represents a plasma or neutral species
/// Can be time evolving, but the default base class
/// is not evolving.
class Species {
public:
  Species() {}
  virtual ~Species() {}
  
  BoutReal Nnorm, Tnorm, Vnorm; /// Normalisation factors
  
  Field3D N;                    // Density
  Field3D P;                    // Pressure
  Field3D NV;                   // Momentum
  
  Field3D T;                    // Temperature
  Field3D V;                    // Velocity
  
  /// Calculate time derivatives
  /// overwriting any previous values
  virtual void evolve(BoutReal UNUSED(time)) {
    // Zero time derivatives, since these may be used for reactions
    ddt(N) = 0.0;
    ddt(P) = 0.0;
    ddt(NV) = 0.0;
  }
  
};

/// Map from string to Species pointers
using SpeciesMap = std::map<std::string, Species*>;

class FluidSpecies : public Species {
public:
  FluidSpecies(std::string name,
               Options *opt,
               Solver *solver,
               bool restarting);
  
  void evolve(BoutReal UNUSED(time)) override;
  
private:
  std::string name;
  
  bool bndry_flux_fix;
  BoutReal gamma_sound; // Ratio of specific heats in numerical dissipation term
  
  int density_sheath;  // How to handle density boundary?
  int pressure_sheath; // How to handle pressure boundary?
  BoutReal sheath_gamma; // Sheath heat transmission factor

  BoutReal anomalous_D, anomalous_chi; // Anomalous transport

  BoutReal viscos; // Numerical dissipation terms
  
  Field3D eta_i;        // Braginskii ion viscosity
  bool ion_viscosity;   // Braginskii ion viscosity on/off

  BoutReal Nnorm, Tnorm, Omega_ci, Cs0;
  BoutReal tau_e0;
};

#endif // __SPECIES_HXX__
