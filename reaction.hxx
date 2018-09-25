/// Defines interface to volume reaction rates

#ifndef __REACTION_HXX__
#define __REACTION_HXX__

#include "bout/generic_factory.hxx"
#include "field3d.hxx"
#include "unused.hxx"

#include <map>
#include <string>

#include "species.hxx"

/// Map from string (species name) to volume source
using SourceMap = std::map<std::string, Field3D>;

/// Represent a reaction
class Reaction {
public:
  virtual ~Reaction() {}
  
  /// Update all species properties
  ///
  /// @param[in] Tnorm   Temperature [eV]
  /// @param[in] Nnorm   Density [m^-3]
  /// @param[in] Vnorm   Velocity [m/s]
  /// @param[in] Freq    Frequency [s^-1]
  virtual void updateSpecies(const SpeciesMap &UNUSED(species),
                             BoutReal UNUSED(Tnorm), BoutReal UNUSED(Nnorm),
                             BoutReal UNUSED(Vnorm), BoutReal UNUSED(Freq)) {}

  /// Return density source for a set of species
  virtual SourceMap densitySources() { return {}; }

  /// Return momentum sources
  virtual SourceMap momentumSources() { return {}; }

  /// energy sources
  virtual SourceMap energySources() { return {}; }

  /// Return a description
  virtual std::string str() const { return "Unknown reaction"; }
};

/// Interpolation of quantities onto cell edges
struct InterpCell {
  InterpCell(const Field3D &f, const DataIterator &i) {
    c = f[i];
    l = 0.5 * (f[i.ym()] + c);
    r = 0.5 * (c + f[i.yp()]);
  }
  InterpCell(const Field2D &f, const DataIterator &i) {
    c = f[i];
    l = 0.5 * (f[i.ym()] + c);
    r = 0.5 * (c + f[i.yp()]);
  }

  BoutReal l, c, r;
};


#endif // __REACTION_HXX__
