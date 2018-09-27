/// Defines interface to volume reaction rates

#ifndef __REACTION_HXX__
#define __REACTION_HXX__

#include "bout/generic_factory.hxx"
#include "field3d.hxx"
#include "unused.hxx"
#include "options.hxx"

#include <map>
#include <string>

#include "species.hxx"
#include "bout/coordinates.hxx"
#include "macro_for_each.hxx"

/// Map from string (species name) to volume source
using SourceMap = std::map<std::string, Field3D>;

/// Represent a reaction
class Reaction {
public:
  Reaction() {}
  Reaction(Options *UNUSED(options)) {}
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

using ReactionFactory = Factory<Reaction, std::function<Reaction *(Options*)>>;

// We need to specialise the helper class to pass arguments to the constructor
template<typename DerivedType>
class RegisterInFactory<Reaction, DerivedType> {
public:
  RegisterInFactory(const std::string &name) {
    ReactionFactory::getInstance().add(
        name, [](Options *option) -> Reaction * { return new DerivedType(option); });
  }
};

//////////////////////////////////////////////////////////////////////////////
// Tools for integrating reaction rates over cells

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
  
  std::array<BoutReal,3> weights; // c l r

  int size() { return weights.size(); }
  BoutReal operator[](int index) const { return weights[index]; }
};

#define _QUADVALUE_DEFINE(var) \
  for(QuadValue _qv_ ## var (var, _int_i); !_control_outer; )

#define _BOUTREAL_DEFINE(var) \
  for(BoutReal var = _qv_ ## var [_int_j]; !_control_inner; )

/// Loop over cells in a region, and a set of quadrature points
/// in those cells, to integrate one or more expressions.
/// These weights sum to 1, resulting in volume averaged values.
///
/// Example:
///
/// Field3D f = ..., g = ...; // Input fields
///
/// Field3D result = 0.0;
/// CELL_AVERAGE(i, result.region(RGN_NOBNDRY), // cell index
///              mesh->coordinates(),   // Coordinates pointer
///              weight,   // Quadrature weight variable 
///              f, g) {   // Fields to interpolate to quadrature points
///
///   // This will be evaluated several times for each cell index i
///   result[i] += weight * ( f + g ); // f and g here are BoutReal
/// }
/// 
/// Introduces the following local variables which should not be used:
///   _int_i, _int_j                  Index and loop counters
///   _control_outer, _control_inner  Booleans to control loops
///   _QR                             QuadRule object with weights
/// In addition, for each input field a new variable with
/// "_" prepended will be defined for the QuadValue object.
///
#define CELL_AVERAGE(indx, region, coord, weight, ...)                  \
  /* Iterate over region with internal index */                         \
  for (auto _int_i : region)                                            \
    /* Introduce outer loop control variable _control_outer */          \
    /* The if(false);else pattern is used to avoid          */          \
    /* accidentally capturing an 'else' after the macro     */          \
    if(bool _control_outer = false) ; else                              \
      for (auto indx = _int_i; !_control_outer;)                        \
        for ( QuadRule _QR(coord, _int_i); !_control_outer;)            \
          /* Define QuadValue variables, interpolating */               \
          MACRO_FOR_EACH(_QUADVALUE_DEFINE, __VA_ARGS__)                \
            /* Loop over weights */                                     \
            for (int _int_j = 0; !_control_outer; )                     \
              /* Setting _control_outer=true causes enclosing loops to exit */ \
              for (_control_outer = true; _int_j < _QR.size(); ++_int_j) \
                /* Define BoutReal variables matching original var names */ \
                if (bool _control_inner = false) ; else                 \
                  for ( BoutReal weight = _QR[_int_j]; !_control_inner; ) \
                    MACRO_FOR_EACH(_BOUTREAL_DEFINE, __VA_ARGS__)       \
                      /* Innermost expression. Run once, setting */     \
                      /* _control_inner=true so enclosing loops exit */ \
                      /* Note: The single '=' is intentional here */    \
                      if (!(_control_inner = true)) ; else              \


#endif // __REACTION_HXX__
