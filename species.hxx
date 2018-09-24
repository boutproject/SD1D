#ifndef __SPECIES_HXX__
#define __SPECIES_HXX__

#include "bout/generic_factory.hxx"

#include <map>
#include <string>

class Species {
public:
  BoutReal Nnorm, Tnorm, Vnorm; /// Normalisation factors
  Field3D N;  // Density
  Field3D P;  // Pressure
  Field3D NV; // Momentum
  Field3D T;  // Temperature
  Field3D V;  // Velocity
};

/// Map from string to Species
using SpeciesMap = std::map<std::string, Species>;

#endif // __SPECIES_HXX__
