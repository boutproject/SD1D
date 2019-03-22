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
#include <cmath>

/// Represents a plasma or neutral species
/// Can be time evolving, but the default base class
/// is not evolving.
class Species {
public:
  Species() : AA(nan("")), ZZ(nan("")) {}
  Species(BoutReal AA, BoutReal ZZ) : AA(AA), ZZ(ZZ) {}
  virtual ~Species() {}

  BoutReal AA; // Atomic mass e.g. Deuterium = 2
  BoutReal ZZ; // Charge e.g. Hydrogen = +1, electron = -1
  
  BoutReal Nnorm, Tnorm, Vnorm; /// Normalisation factors
  
  Field3D N;                    // Density
  Field3D P;                    // Pressure
  Field3D NV;                   // Momentum, includes factor of AA
  
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
  FluidSpecies(std::string name, Options *opt, Solver *solver,
               Datafile &restart, bool restarting, BoutReal Nnorm,
               BoutReal Tnorm, BoutReal Omega_ci, BoutReal Cs0);

  void evolve(BoutReal time) override;
  
private:
  std::string name;
  
  bool bndry_flux_fix;
  BoutReal gamma_sound; // Ratio of specific heats in numerical dissipation term

  bool sheath_outflow; // True if species flows out at >= sound speed
  int density_sheath;  // How to handle density boundary?
  int pressure_sheath; // How to handle pressure boundary?
  BoutReal sheath_gamma; // Sheath heat transmission factor

  BoutReal anomalous_D, anomalous_chi; // Anomalous transport

  BoutReal viscos; // Numerical dissipation terms
  
  Field3D eta_i;        // Braginskii ion viscosity
  bool ion_viscosity;   // Braginskii ion viscosity on/off

  BoutReal Nnorm, Tnorm, Omega_ci, Cs0;
  BoutReal tau_e0;
  
  Field2D NeSource, PeSource; // Volume sources
  Field2D NeSource0;          // Used in feedback control

  // Upstream density controller
  BoutReal density_upstream; // The desired density at the lower Y (upstream)
                             // boundary
  BoutReal density_controller_p, density_controller_i; // Controller settings
  bool density_integral_positive; // Limit the i term to be positive
  bool density_source_positive;   // Limit the source to be positive

  BoutReal density_error_lasttime,
      density_error_last;          // Value and time of last error
  BoutReal density_error_integral; // Integral of error
};

#endif // __SPECIES_HXX__
