/// Defines interface to volume reaction rates

#ifndef __REACTION_HXX__
#define __REACTION_HXX__

#include "bout/generic_factory.hxx"

#include <map>
#include <string>

#include "species.hxx"

/// Represent a reaction 
class Reaction {
public:
  Field3D &rate(const SpeciesMap &species) {
    if (!cache_valid) {
      cache = calculate(species);
      cache_valid = true;
    }
    return cache;
  }
  
protected:
  /// Calculate the rate from given species information
  virtual Field3D calculate(const SpeciesMap &species) = 0;

private:
  bool cache_valid = false; ///< Is the cached value valid?
  Field3D cache; ///< The calculated value
};

/// Map from string to Reaction
using ReactionMap = std::map<std::string, Reaction>;

#endif // __REACTION_HXX__
