/// Defines interface to volume reaction rates

#ifndef __REACTION_HXX__
#define __REACTION_HXX__

#include "bout/generic_factory.hxx"
#include "field3d.hxx"

#include <map>
#include <string>

#include "species.hxx"

/// Map from string (species name) to volume source
using SourceMap = std::map<std::string, Field3D>;

/// Represent a reaction
class Reaction {
public:
  /// Update all species properties
  virtual void updateSpecies(const SpeciesMap &species, BoutReal Tnorm,
                             BoutReal Nnorm, BoutReal Vnorm) {}

  /// Return density source for a set of species
  virtual SourceMap densitySources() { return {}; }

  /// Return momentum sources
  virtual SourceMap momentumSources() { return {}; }

  /// energy sources
  virtual SourceMap energySources() { return {}; }
};

/// Map from string to Reaction
using ReactionMap = std::map<std::string, Reaction*>;

#endif // __REACTION_HXX__
