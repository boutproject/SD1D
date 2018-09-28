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

  /// Initialise time integration
  virtual void init(bool UNUSED(restarting), Solver *UNUSED(solver)) {}
  /// Calculate time derivatives
  /// overwriting any previous values
  virtual void evolve(BoutReal UNUSED(time)) {}
  
};

/// Map from string to Species
using SpeciesMap = std::map<std::string, Species>;

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
  
};

#endif // __SPECIES_HXX__
