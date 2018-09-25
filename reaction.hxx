/// Defines interface to volume reaction rates

#ifndef __REACTION_HXX__
#define __REACTION_HXX__

#include "bout/generic_factory.hxx"
#include "field3d.hxx"
#include "unused.hxx"

#include <map>
#include <string>

#include "species.hxx"
#include "bout/coordinates.hxx"

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

/// Values to be used in a quadrature rule (Simpson's)
struct QuadValue {
  QuadValue(const Field3D &f, const DataIterator &i) {
    values[0] = f[i];
    values[1] = 0.5 * (f[i.ym()] + values[0]);
    values[2] = 0.5 * (values[0] + f[i.yp()]);
  }
  QuadValue(const Field2D &f, const DataIterator &i) {
    values[0] = f[i];
    values[1] = 0.5 * (f[i.ym()] + values[0]);
    values[2] = 0.5 * (values[0] + f[i.yp()]);
  }

  BoutReal values[3];   // c l r

  BoutReal operator[](int index) const {return values[index];}
};

/// Simpson's rule weights
struct QuadRule {
  QuadRule(Coordinates *coord, const DataIterator &i) {
    Field2D J = coord->J;
    auto ym = i.ym();
    auto yp = i.yp();
    
    weights[0] = 4./6.;
    weights[1] = 0.5*(J[ym] + J[i]) / (6. * J[i]);
    weights[2] = 0.5*(J[yp] + J[i]) / (6. * J[i]);
  }
  
  BoutReal weights[3]; // c l r

  BoutReal operator[](int index) const {return weights[index];}

  int indices[3] = {0,1,2};
};

#endif // __REACTION_HXX__
